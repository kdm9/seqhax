TOOLS=(
    'pecheck'
    )

cd $(dirname $(readlink -f $0))
for subtool in "${TOOLS}"
do
    cd $subtool
    bash test.sh
    cd ..
done

