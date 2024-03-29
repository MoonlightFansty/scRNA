# Method Comparation
## 一、归一化
由于单细胞测序每个步骤中固有的可变性，相同细胞的计数深度可能会有所不同。因此，当基于计数数据比较细胞之间的基因表达时，任何差异都可能仅由于采样效应而产生。标准化可以通过缩放计数数据来获得细胞间正确的相对基因表达丰度来解决这个问题。这里我们对使用最多的两种方法：Seurat标准scale和SCTransform进行比较。提取出部分单细胞表达数据，对每个基因做简单的统计包括基因表达值均值、方差以及可检测率。A中可以看到，在均值较小的一定范围内，均值和方差有线性关系。均值再变大时，方差呈现过度分散。B是基因表达均值和基因可检测率的关系图，图中红线是假设基因表达均值符合泊松分布的期望可检测率。可以看到高表达和低表达的基因符合预期值，基因表达量局中的基因的可检测率低于预期。 \
![](https://github.com/MoonlightFansty/scRNA/blob/main/Method%20Comparation/Figures/1-1.png)
在图A中，使用SCTransform方法标准化，并看一下基因表达量的改变。我们挑选了三个基因MALAT1，RPL10和FTL。上半部显示矫正前，下半部显示矫正后基因表达量和测序深度的关系。可以看到经过SCTransform我们达到了目的，基因表达量不再和每个细胞的测序深度有相关性。
同时使用Seurat标准scale方法，为保证可比性，我们regress out总基因表达量，如图B。从标准流程出来的残差和细胞总表达量的关系中，可以看到基因表达量同样的不再和每个细胞的测序深度有相关性，但数据离散性很大。而对于像FTL这样的基因，用SCTransform后能明显看到其表达分为高和低两组细胞，而用标准流程此效果不佳，因此选择SCTransform方法。 \
![](https://github.com/MoonlightFansty/scRNA/blob/main/Method%20Comparation/Figures/1-2.png)

## 二、去除批次效应
接着，我们比较了多种批次整合方法并评估了它们带来的影响。经典的数据整合方法包括CCA，MNN、Scanorama、scGen、LIGER、BBKNN和Harmony等。这里我们参考了2021年的一篇单细胞批次处理的基准测评综述。从图中可知，在标准化方法中SCTransform性能较好，这与我们验证的结果相同，同时在整合方法中harmony性能占优。 \
![](https://github.com/MoonlightFansty/scRNA/blob/main/Method%20Comparation/Figures/2-1.png)
由于新方法的出现，我们使用了LISI评分，比较harmony与其他几种方法整合效果。发现在现有工具中harmony的整合性能依然优越，并且运行速度很快。 \
![](https://github.com/MoonlightFansty/scRNA/blob/main/Method%20Comparation/Figures/2-2.png) \
最后，对harmony算法的应用进行评估，观察其整合效果，并考量生物学意义是否丢失。首先，部分样品的单细胞转录组数据整合后，如果不使用harmony等算法去除样品差异，默认的降维聚类分群。图A，我们对细胞类型进行注释，并绘制UMAP图。在运行harmony后进行聚类，同样进行UMAP展示，如图B所示。我们发现因为harmony算法，上皮细胞被打散到了其它免疫细胞里面，说明出现了过度校正，导致生物学差异丢失。 \
![](https://github.com/MoonlightFansty/scRNA/blob/main/Method%20Comparation/Figures/2-3.png)
