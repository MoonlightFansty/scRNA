# Seurat pipeline
## Preprocessing & Cluster
我们使用Seurat包(版本4.1.1)进行后续分析，读入稀疏矩阵并转换为Seurat对象，合并所有单细胞基因表达数据。接着，由于低质量细胞或空液滴通常只能检测到非常少的基因，两个或多个细胞被同时捕获通常会有很高的基因数，低质量和濒死细胞常表现出广泛的线粒体污染，所以需要过滤转录组数据。根据图1-1A、B确定阈值，筛选出基因数为500-10，000，UMI计数为1000-100，000、线粒体读数的分数小于30%的细胞。根据图1-1C、D可以看到经过滤后，随着测序深度的增加，测得特征基因和UMI计数相关性更强，说明去除了低质量的细胞、空液滴和含多个细胞的液滴。线粒体基因比例低且分布均匀，说明过滤掉了低质量、濒死以及破损的细胞。 \
![1-1](https://github.com/MoonlightFansty/scRNA/blob/main/Seurat/Figures/1-1.png) \
过滤后，使用具有3000个可变特征的SCTransform对UMI计数进行方差稳定，同时回归出UMI的总数和线粒体转录物含量的比例。为了降低单细胞数据集的维数，在集成数据矩阵上进行了主成分分析（PCA）。选取最优的主成分数需要满足主成分累积贡献大于90%，PC本身对方差贡献小于5%，两个连续PCs之间差异小于0.1%。图1-2A在Eblow图中，经13个PC之后，显著性大幅下降，也就是前13个维度包含了大部分的样本信息。如图1-2B所示，在Heatmap图中，以PCA分数排序的基因和细胞可视化选定PC的最变异基因。观察不同PC的重要性，查看基因差异是否明显，找到热图开始区分不明显的PC数。综合以上方法，选择15个主成成分作为参数用于后续分析。 \
![1-2](https://github.com/MoonlightFansty/scRNA/blob/main/Seurat/Figures/1-2.png) \
主细胞簇通过Seurat提供的FindClusters功能进行识别，分辨率的设置由clustrees图和ROGUE图确定。如图1-3A所示，通过这个图我们会发现当resolution大于0.2后，不同cluster在不同的分辨率下会存在越来越多的相互交织，这暗示着我们可能存在分群过度的情况了，所以在这里我们可以选择0.2作为初步的分群resolution来进行后面的分析。图1-3B使用ROUGE验证我们的分群效果，一个细胞类群的ROGUE值表征了该细胞类群的“纯度”：越接近于1表示细胞类群越“纯”。可以看到分群效果很好，因此确定resolution为0.2。 \
![1-3](https://github.com/MoonlightFansty/scRNA/blob/main/Seurat/Figures/1-3.png) \
然后用2D的tSNE或UMAP图将它们可视化。tSNE图的细胞聚类更为紧凑，但不能反映全剧结构。UMAP图中除聚类外，整体结构反映了细胞间的真实距离，计算速度更快。由于细胞数量较多，使用UMAP图进行结果展示。利用 FindMarkers 命令，可以找到找到各个细胞类型中与其他类别的差异表达基因，作为该细胞类型的生物学标记基因。
![1-4](https://github.com/MoonlightFansty/scRNA/blob/main/Seurat/Figures/1-4.png) \
接着，根据先前的研究，获取了肺腺癌主要分类的生物标记基因。通过DimPlot图对19个细胞类群进行注释，得到单细胞三个主要类群（上皮、免疫、基质细胞）。并查看每个样本中，三类细胞所占的比例，上皮细胞在肿瘤患者中占比更高，符合生物学常识。 \
![1-5](https://github.com/MoonlightFansty/scRNA/blob/main/Seurat/Figures/1-5.png) \
然后进行细胞周期影响的判断，当这群细胞占比很少时，可以不做处理；占比较大，或者处于不同细胞周期时间点(G1, S, G2, M)的细胞相互分开，这时就需要消除细胞周期的影响，留下有意义的生物学差异。从图1-6A中，我们可以看出S.Score较高的为S期，G2M.Score较高的为G2M期，都比较低的为G1期。而从图1-6B可以看到细胞周期不同时间点的细胞混合得很均匀，对聚类结果的影响很小。 \
![1-6](https://github.com/MoonlightFansty/scRNA/blob/main/Seurat/Figures/1-6.png) \
最后，我们查看了组织类型和患者的细胞分布状况，不同组织类型或患者的单细胞转录组混合在许多簇中，不包括一般的批次效应。部分形成的肿瘤或患者特异性簇，表明存在的生物学差异。 \
![1-7](https://github.com/MoonlightFansty/scRNA/blob/main/Seurat/Figures/1-7.png)

## 二、Recluster
在获得了三个主群分类后，我们在接着对每个主类群进行亚群重聚类。细胞类型的基因marker改编自Habermann以及Tata and Rajagopal等人。Vieira Braga和Travaglini等人的细胞类型marker用于验证手动细胞类型注释。
首先，对上皮细胞进行重聚类，绘制所有上皮单细胞转录组的前20个主成分的UMAP，鉴定出了基底上皮细胞、肺泡1型细胞、肺泡2型细胞、细支气管外分泌细胞、纤毛上皮细胞和神经内分泌细胞等不同簇。 \
![2-1](https://github.com/MoonlightFansty/scRNA/blob/main/Seurat/Figures/2-1.png) \
通过标记基因和基因特征鉴定了内皮和淋巴内皮细胞，成纤维细胞，肌成纤维细胞和平滑肌细胞以及间皮细胞的不同簇。 \
![2-2](https://github.com/MoonlightFansty/scRNA/blob/main/Seurat/Figures/2-2.png) \
最后，我们根据典型的标记基因和基因特征鉴定了组织驻留和单核细胞来源的巨噬细胞，单核细胞，髓细胞和浆细胞样树突状细胞，肥大细胞以及T，NK，B和浆细胞的不同簇。
![2-3](https://github.com/MoonlightFansty/scRNA/blob/main/Seurat/Figures/2-3.png)
