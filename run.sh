# $1 is w
# $2 is initial guess of shooting parameter
# $3 is r1
# $4 is rmax

echo 1.0 > parameters
echo 1.0 >> parameters
echo 0.9 >> parameters
echo $1 >> parameters
echo $2 >> parameters
echo $3 >> parameters
echo $4 >> parameters

python run.py

mkdir -p data/w=$1
mv {pars,integral-data,functions}.dat data/w=$1/
mv {phi,psi}.pdf data/w=$1/
mv parameters data/w=$1/
