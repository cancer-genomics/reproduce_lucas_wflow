Downloading TCGA cn data via rtcga

library(RTCGA)
releaseDate="2016-01-28"
downloadTCGA(cancerTypes = "LUSC", destDir = ".", date = releaseDate, dataSet = "Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3" ), silent=TRUE)
downloadTCGA(cancerTypes = "LUAD", destDir = ".", date = releaseDate, dataSet = "Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3" ), silent=TRUE)


# 7cd84ea888f4eb5a32ee72c44aef66c6 LUSC.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt
# e5b29b5b5e335816f103f3158c2b7755 LUAD.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt
