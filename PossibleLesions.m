function ind=PossibleLesions(minn,maxx)

load('FeatureFile.mat');

ind1=(a>=minn(1) & a<=maxx(1)); 
ind2=(per>=minn(2) & per<=maxx(2));
ind3=(MrM>=minn(3) & MrM<=maxx(3));
ind4=(circularity>=minn(4) & circularity<=maxx(4));
ind5=(i_green>=minn(5) & i_green<=maxx(5));
ind6=(i_sc>=minn(6) & i_sc<=maxx(6));
ind7=(m_green>=minn(7) & m_green<=maxx(7));
ind8=(m_sc>=minn(8) & m_sc<=maxx(8));
ind9=(NI_green>=minn(9) & NI_green<=maxx(9));
ind10=(NI_sc>=minn(10) & NI_sc<=maxx(10));
ind11=(NM_green>=minn(11) & NM_green<=maxx(11));
ind12=(NM_sc>=minn(12) & NM_sc<=maxx(12));
ind13=(IDarkest>=minn(13) & IDarkest<=maxx(13));
ind14=(FiltSig1>=minn(14) & FiltSig1<=maxx(14));
ind15=(FiltSig2>=minn(15) & FiltSig2<=maxx(15));
ind16=(FiltSig4>=minn(16) & FiltSig4<=maxx(16));
ind17=(FiltSig8>=minn(17) & FiltSig8<=maxx(17));
ind18=(StdFiltSig1>=minn(18) & StdFiltSig1<=maxx(18));
ind19=(StdFiltSig2>=minn(19) & StdFiltSig2<=maxx(19));
ind20=(StdFiltSig4>=minn(20) & StdFiltSig4<=maxx(20));
ind21=(StdFiltSig8>=minn(21) & StdFiltSig8<=maxx(21));
ind22=(MaxCorr>=minn(22) & MaxCorr<=maxx(22));
ind23=(MinCorr>=minn(23) & MinCorr<=maxx(23));
ind24=(AvgCorr>=minn(24) & AvgCorr<=maxx(24));
ind26=(diffr>=minn(26) & diffr<=maxx(26));
ind27=(diffg>=minn(27) & diffg<=maxx(27));
ind28=(diffb>=minn(28) & diffb<=maxx(28));
ind29=(diffh>=minn(29) & diffh<=maxx(29));
ind30=(MajorAxis>=minn(30) & MajorAxis<=maxx(30));
ind31=(MinorAxis>=minn(31) & MinorAxis<=maxx(31));

indd=ind1&ind2&ind3&ind4&ind5&ind6&ind7&ind8&ind9&ind10&ind11&ind12&ind13&ind14&ind15&ind16&ind17;
ind=indd&ind18&ind19&ind20&ind21&ind22&ind23&ind24&ind26&ind27&ind28&ind29&ind30&ind31;