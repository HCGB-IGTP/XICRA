## get path where is conda environment
path=`which XICRA`
dir=`dirname $path`

mkdir tmp
cd tmp

## get sRNAtoolboxDB
echo ""
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++ "
echo "## Downloading sRNAtoolboxDB"
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++ "
wget http://bioinfo2.ugr.es/sRNAtoolboxDB/sRNAtoolboxDB.tgz
tar -zxvf sRNAtoolboxDB.tgz
echo ""

echo ""
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++ "
echo "## Downloading sRNAbench.jar"
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++ "
wget http://bioinfo2.ugr.es/sRNAtoolboxDB/exec/sRNAbench.jar
mv sRNAbench.jar sRNAtoolboxDB/exec/
echo ""

echo ""
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++ "
echo "## Downloading sRNAde.jar "
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++ "
wget http://bioinfo2.ugr.es/sRNAtoolboxDB/exec/sRNAde.jar
mv sRNAde.jar sRNAtoolboxDB/exec/
echo ""

echo ""
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++ "
echo "## Downloading miraligner.jar"
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++ "
wget https://github.com/lpantano/seqbuster/raw/miraligner/modules/miraligner/miraligner.jar
echo ""


echo ""
## mv files to conda enviornment folder
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++ "
echo "Moving files"
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++ "
ln -s sRNAtoolboxDB/exec/sRNAbench.jar
mv sRNAtoolboxDB $dir
mv sRNAbench.jar $dir
mv miraligner.jar $dir

echo ""
echo ""
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++ "
echo "Remove tmp files"
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++ "
cd ../
rm -r tmp

echo ""
echo ""
echo "Done..."
echo ""

