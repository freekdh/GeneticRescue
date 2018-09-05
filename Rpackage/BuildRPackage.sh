rm -rf pkgIntrogression
R CMD BATCH --vanilla --no-timing MakePackage.R /dev/null

# Annoying thing with Mac and Linux. Mac needs the ''.
sed -i '$d' pkgIntrogression/DESCRIPTION
sed -i '$d' pkgIntrogression/DESCRIPTION
echo "Imports: Rcpp (>= 0.12.17), BH, RcppProgress" >> pkgIntrogression/DESCRIPTION
echo "LinkingTo: Rcpp, BH, RcppProgress" >> pkgIntrogression/DESCRIPTION

# move shiny files to inst folder
mkdir -p pkgIntrogression/inst/ShinyApp
cp ./ShinyApp/server.r ./pkgIntrogression/inst/ShinyApp
cp ./ShinyApp/ui.r ./pkgIntrogression/inst/ShinyApp

R CMD INSTALL pkgIntrogression