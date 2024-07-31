## Update 20240731

### CrossTalkeR install

まずはCellPhoneDBのinstall(リポジトリ内に最新版v5.00あり)

```bash
conda create -n cpdb python=3.8
conda activate cpdb

pip install cellphonedb
cellphondbは手動でダウンロード
```

<使用手順>  
１：mail.shの以下を実行する環境に合わせて修正すること  
  cellphonedb_proc.R, shaping_hcount_meta.py, cellphonedb_res_shaping.py.   
  のそれぞれのファイルディレクトリの指定部分.   

2 ：shaping_hcount_meta.pyの以下を実行する環境に合わせて修正すること.   
   cellphonedb.zipのディレクトリの指定  
   
3 ：main.shを用いたコマンドを実行  
<例>  
bash main.sh \  
-f file_path \ ファイルパスの指定  
-p project_name \　プロジェクト名等を指定  
-o output_dir \ 結果の出力先(パス)を指定  
-t statistical_analysis_significant_means_Text_file_name \ statistical_analysis_significant_means_Textの名前を指定した名前にリネーム  
-a mail_address 実行が完了したことを知らせるメールアドレス先  
