#RおよびBioconductorを使ったバイオデータ解析

講習では【実習】を皆さんとやっていきます。早くできてしまった人は、【発展】をやってみたりしてください。Rのインストールは以下の統合TVを参考に行っておいてください。

* 参考:[統合TV](http://togotv.dbcls.jp/) 
![togotv](http://togotv.dbcls.jp/images/togotv_banner.png)
 * [統計解析ソフト「R」の使い方 導入編](http://togotv.dbcls.jp/20090313.html)
 * [統計解析ソフト「R」の使い方 導入編(MacOSX)](http://togotv.dbcls.jp/20090618.html)
 * [統計解析ソフト「R」の使い方 正規化編](http://togotv.dbcls.jp/20090319.html)
 * [統計解析ソフト「R」の使い方 ヒートマップ編](http://togotv.dbcls.jp/20091219.html)
 * [統計解析ソフト「R」での立廻り](http://togotv.dbcls.jp/20111107.html)

--
## はじめに

今回使うRのコードやデータはgithubのサイトにあります。そこのデータをまとめてダウンロードするには、UNIXのコマンドラインで

    git clone https://github.com/bonohu/AJACS52
を実行すれば現在居るディレクトリ(current working directory)にAJACS52というディレクトリが作成され、それ以下にgithubにあるファイルがすべてダウンロードされますが、今回の講習会の環境ではできませんので、[https://github.com/bonohu/AJACS52](https://github.com/bonohu/AJACS52)の画面右下の方の「Download ZIP」のリンクをたどると一括ダウンロードできます。

Rの実行に関しては主に2種類あって
 * R をコマンドラインから使う場合
 * Rのアイコンをダブルクリックして使う場合
  * [Rstudio](http://www.rstudio.com/)からRを実行する 

があります。いずれの使用においても、`setwd()`コマンドを使って今日の講習で使うディレクトリ名を指定します。
**これが全員できるまでここから先には進みません。**

## 1. pie chart(パイチャート)
【実習1】Rを起動し以下のファイルに記述されたRのコマンド(`01piechart.r`)を実行しなさい。データファイルとして必要な`srabystudy.txt`をgithubのサイトからダウンロードし、現在作業しているディレクトリにおいておく必要があります。実行後、そのディレクトリに`pie1.png`というファイルが新たに生成されていることを確認しなさい。

    png("pie1.png")
    dat <- read.delim2("srabystudy.txt", header=F)
    names(dat) <- c("Study", "freq")
    dat <- cbind(dat, serial=seq(dim(dat)[1]))
    dat2 <- dat[ dat$freq != 0, ]
    pie(dat2$freq, labels=paste(dat2$Study, dat2$freq), col=rainbow(dim(dat[1])[dat2$serial]), main="SRA by Study")
    dev.off()

先ほど説明のあったSRA(Sequence Read Archive)のメタデータを取りまとめてStudyのタイプで分けた結果がこのデータ(`srabystudy.txt`)で、それを元にpie chartを描いてみるのが本課題でした。

【発展1】上記の`srabystudy.txt`を別のデータに変えて実行してみましょう。
 * `srabystudy_orig.txt` : SRAのStudyの分類をoriginalデータ
 * `srabyorganism.txt` : SRAに登録されている生物種で分類したデータ

> これらの結果は、お配りしたファイル群の`results`というフォルダの中においてあります(それぞれ、`pie1.png`,`pie2.png`,`pie3.png`)。実は、このスクリプトは[Bono H et al. Genome Res 2003](http://genome.cshlp.org/content/13/6b/1318.full)の図1を作成するために使ったものです。実際には、このRコードを生成するPerlのプログラムを作成してからそれをRで実行していました。現在よく使われているGSEA(Gene Set Enrichment Analysis)のハシリで、遺伝子にアノテーションされたGene Ontologyの高次(根っこに近い)タームで集計してその分布を見ていました

## 2. hierachical clustering(階層的クラスタリング)
【実習2】Rを起動し以下のファイルに記述されたRのコマンド(`02hclust.r`)を実行しなさい。データファイルとして必要な`matrix.txt`をgithubのサイトからダウンロードし、現在作業しているディレクトリにおいておく必要があります。実行後、そのディレクトリに`hclust.png`というファイルが新たに生成されていることを確認しなさい。

    png("hclust.png")
    d <- read.table('matrix.txt')
    c <- hclust(as.dist(d), method = 'average')
    plot(c, hang=-1)
    dev.off()
    
このプログラムはかつての蛋白質・核酸・酵素に寄稿したレビューで紹介したもので、当時階層的クラスタリングをフリーウェアで実行する手段として重宝された(はず)。[Bono H and Nakao MC, PNE 2003](http://www.ncbi.nlm.nih.gov/pubmed/12638180)

## 3. Affymetrixデータ正規化 → 主成分分析(PCA)
### Affymetrixデータ正規化
このパートは[ぼうのブログ: justRMAでnormalize](http://bonohu.jp/blog/2013/06/17/justrma/ "ぼうのブログ: justRMAでnormalize")を参考に。RMAの計算が重く、パソコンによってはメモリ不足となり実行不可能となることが予想されます。

【発展2】Bioconductorを利用しましょう。Rを起動して以下のコマンドでaffy libraryをインストールしなさい。

    source("http://bioconductor.org/biocLite.R")
    biocLite("affy")
    
そして、affy library中のjustRMA関数を利用して、RMA(Robust Multichip Average)正規化してみましょう。自らのAffymetrix Genechip データ(CELファイル)がある人はそれを、ない人はGSE17264 Comparative transcriptome analysis of dedifferentiation in porcine mature adipocytes and follicular granulosa cellsにある生データ(CELファイル)で実行しましょう。実行する際、Session → Set Working Directory を、CELファイルをダウンロードしてきたディレクトリに指定することがポイントです。それができたら、以下のコード(`03justrma.r`)を実行します。

    source("http://bioconductor.org/biocLite.R")
    biocLite("affy")
    library(affy)
    write.exprs(justRMA(), file="RMA.txt")
    
その結果生成される`RMA.txt`ファイルがRMAによって正規化された遺伝子発現データになります。

【発展3】上記の実習のサンプルデータは **GPL3533 [Porcine] Affymetrix Porcine Genome Array** と呼ばれる(カタログ)マイクロアレイ(プラットフォーム)を使ったデータでした。同じマイクロアレイ(プラットフォーム)のデータであればjustRMAを実行することが出来ます。同じプラットフォームのデータの中からiPS細胞のものを探しだして上記のデータに混ぜてjustRMAを実行しなさい。

> 解答例: [GSE15472 Induced Pluripotent Stem Cells from the Pig Somatic Cells](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15472) の3つのサンプル(GSM388154,GSM388155,GSM388156)がそれです。 CELファイルをダウンロードして同じdirectoryに置き、再度justRMAを実行してみましょう

> このマイクロアレイデータは以下の論文で我々が発表したもので、全く同じ手段で正規化を行いました。また、発展課題にある公共遺伝子発現データを足したデータ解釈はやはりこの論文で実際に発表したもので、我々のマイクロアレイデータの生物学的解釈(DFATがMA(Mature AdipocyteよりもiPSに近い)に役だっています。 [Ono H, Oki Y, Bono H, Kano K. BBRC 2011](http://www.ncbi.nlm.nih.gov/pubmed/21419102)

### 主成分分析(PCA)
PCA(Primary Component Analysis)。 このパートは[ぼうのブログ「RでPCA」](http://bonohu.jp/blog/2013/07/13/pcaonr/)を参考に。

【実習3】Rを起動し以下のファイルに記述されたRのコマンド(`03pca.r`)を実行しなさい。データファイルとして必要な`RMA.txt`は【発展2】で作成したものを利用すればよいのですが、現在作業しているディレクトリにおいておく必要があります。  
この発展課題をやらなくても実行できるように用意してあり、`RMA.txt.zip`というファイルをダブルクリックすることで解凍(圧縮を解く)すると`RMA.txt`というファイルが生成されるようになっています。  
PCAを実行するには、以下のRスクリプトを実行します。

    data <- read.table("RMA.txt", header=TRUE, row.names=1, sep="\t", quote="") 
    data.pca <- prcomp(t(data))
    names(data.pca)
    plot(data.pca$sdev, type="h", main="PCA s.d.")
    data.pca.sample <- t(data) %*% data.pca$rotation[,1:2]
    plot(data.pca.sample, main="PCA")  
    text(data.pca.sample, colnames(data), col = c(rep("red", 3), rep("blue",3),rep("green",3),rep("black",3)))
    
図に.CEL.gzの文字がいっぱい入っていて見づらい?そう思った方は[こちらを参照](http://bonohu.jp/blog/2014/09/10/afterjustrma/)

## 4. NGSデータ解析(cuffdiffの結果可視化)
cufflinksパッケージに`cuffdiff`という発現データの差分を計算するプログラムがあります。それを実行したあとのデータを可視化する手段として使うBioconductorのパッケージに`cummeRbund`があります([cuffdiffの使用例](http://dx.doi.org/10.1371%2Fjournal.pone.0104283))。

* 参考: [bonohuの発表資料「NGS解析(RNAseq)」](http://dx.doi.org/10.6084/m9.figshare.1216717)

【発展4】Rを起動し以下のファイルに記述されたRのコマンドを実行しなさい。データファイルとして必要なcuffdiffディレクトリ以下のファイルは手持ちのものか、なければサンプルをUSBメモリでお渡しします。現在作業しているディレクトリにおいておく必要があります。

    source("http://bioconductor.org/biocLite.R")
    biocLite("cummeRbund")
    library("cummeRbund")
    cuff.dir <- "cuffdiff"
    cuff <- readCufflinks(dir=cuff.dir)

準備はここまでで、以下のコマンドで発現密度分布のプロット

    dens <- csDensity(genes(cuff))
    dens

CSV形式でFPKM値を出力などができます。

    gene.matrix <- fpkmMatrix(genes(cuff))
    write.csv(gene.matrix, file="fpkm.csv")
    
## 5.最後に: R Graphical Manualのすゝめ
どんなライブラリがあって、どういった可視化ができるかは[R Graphical Manual](http://rgm3.lab.nig.ac.jp/RGM)が大変便利です。このウェブサイトは国立遺伝学研究所の小笠原理さんによって維持されているRのライブラリに記載されているデモデータをあらかじめ計算して得られるイメージファイルをウェブから閲覧できるというもので、そのライブラリが激しくバージョンアップがなされるRにおいて、現在使用可能かどうかを知る上でも大変有用なものです。

【実習4】ウェブブラウザを起動し、インターネット検索でR Graphical Manualを検索してサイトに辿り着き、サイト内の検索フォームからキーワードcummeRbundで検索しなさい。どのような結果が返ってきたか?前節で紹介したcsDensityのエントリを見つけてみなさい。

【発展5】csDensityのエントリにあるサンプルコードを参考に、csVolcanoやdispersionPlotなどのエントリ中のサンプルコードを現在の例に合うように改変して実行してみなさい。
