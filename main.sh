#!/bin/bash

function usage {
  cat <<EOM
Usage: $(basename "$0") [OPTION]...
  -h          Display help
  -f VALUE    Set a Seurat Object file name f
  -p VALUE    Project Name
  -o VALUE    Output Directory
  -r STRING   (Optional) Set Reshape Directory, if not set, OUTPUT directory will be used
  -t VALUE    Set statistical analysis significant means Text file name t
  -a STRING   Set up an email address to let us know that the analysis is complete.
EOM
}

# 引数別の処理定義
while getopts ":f:p:o:r:t:a:h" optKey; do
  case "$optKey" in
    f) FILENAME="${OPTARG}" ;;
    p) PROJECTNAME="${OPTARG}" ;;
    o) OUTPUT="${OPTARG}" ;;
    r) RESHAPE="${OPTARG}" ;;
    t) TEXTFILE="${OPTARG}" ;;
    a) EMAIL="${OPTARG}" ;;
    h|'-h'|'--help' ) 
      usage 
      exit 0 
      ;;
    * ) 
      usage 
      exit 1 
      ;;
  esac
done

# Conda初期化と環境のアクティブ化
source /home/saito_h/long_read_seq/miniconda3/etc/profile.d/conda.sh
conda activate BioAmadeus

# 変数の表示 (デバッグ用)
echo "FILENAME: $FILENAME"
echo "PROJECTNAME: $PROJECTNAME"
echo "OUTPUT: $OUTPUT"
echo "RESHAPE: ${RESHAPE:-$OUTPUT}"
echo "TEXTFILE: $TEXTFILE"
echo "EMAIL: $EMAIL"

# 出力ディレクトリとリシェイプディレクトリの作成
mkdir -p "$OUTPUT"
mkdir -p "${RESHAPE:-$OUTPUT}"

# Rスクリプトの実行
Rscript "${HOME}/sc_RNAseq/CrosstalkeR/cellphonedb_proc.R" "$FILENAME" "$PROJECTNAME" "$OUTPUT"

# shaping & run cellphoneDB
python "${HOME}/sc_RNAseq/CrosstalkeR/shaping_hcount_meta.py" "${OUTPUT}/${PROJECTNAME}_filtered_hcount.txt" "${OUTPUT}/${PROJECTNAME}_filtered_meta.txt" "$OUTPUT"

# リネーム処理
for file in "${RESHAPE:-$OUTPUT}"/*; do
  if [[ $file == *statistical_analysis_significant_means_* ]]; then
    mv "$file" "${RESHAPE:-$OUTPUT}/$TEXTFILE"
  fi
done

# Reshape
python "${HOME}/sc_RNAseq/CrosstalkeR/cellphonedb_res_shaping.py" "${RESHAPE:-$OUTPUT}" "${RESHAPE:-$OUTPUT}/$TEXTFILE"

# 必要に応じてメール通知の処理を追加
if [ -n "$EMAIL" ]; then
  echo "Analysis complete" | mail -s "Analysis Notification" "$EMAIL"
fi



