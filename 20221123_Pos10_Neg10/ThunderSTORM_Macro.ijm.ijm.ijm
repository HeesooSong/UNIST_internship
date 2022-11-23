open("C:/Users/user/Desktop/20221123_Pos10_Neg10/Negative/PB1484_09/S2_W5_Px5_Py1_C1.tif");
selectWindow("S2_W5_Px5_Py1_C1.tif");
run("Camera setup");
run("Run analysis", "filter=[Wavelet filter (B-Spline)] scale=2.0 order=3 detector=[Local maximum] connectivity=8-neighbourhood threshold=std(Wave.F1) estimator=[PSF: Integrated Gaussian] sigma=1.6 fitradius=3 method=[Weighted Least squares] full_image_fitting=false mfaenabled=false renderer=[Averaged shifted histograms] magnification=5.0 colorizez=false threed=false shifts=2 repaint=50");
selectWindow("Averaged shifted histograms");
close();

run("Export results", "filepath=C:\\Users\\user\\Desktop\\20221123_Pos10_Neg10\\Negative\\PB1484_09\\C1_result.csv fileformat=[CSV (comma separated)] sigma=true intensity=true chi2=true offset=true saveprotocol=false x=true y=true bkgstd=true id=true uncertainty=true frame=true");
close();
