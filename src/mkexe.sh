compile()
{
    code=$1
    allinput="$@"
    echo "================== Compiling "$code" with "$allinput" ================="

    rm -f $code
    #the order of -l has to be tested, it makes a difference in some machine
    for ii in $code.{c,cxx,C,cpp}
    do
        if [ -e $ii ]
        then
            #improtant to use -Wl,--no-as-needed !!
    cmd="g++ $ii  \
        -g -O3 -Wall -Werror \
        $(root-config --auxcflags | tr -s ' ' '\n' | grep c++1)  -I$(root-config --incdir)  \
        -Wl,--no-as-needed -L./ $(root-config --libs) \
        -lTree  -lMinuit -lASImage $lib \
        -o $allinput"
    echo
    echo $cmd
    eval $cmd
    echo
        fi
    done
}

if [ $# == 0 ]
then
    echo no input!!
    exit
fi

compile "$@"

if [ -e $1 ]
then
    echo
    echo checked: $0 $1 done!
    echo

else
    echo
    echo "*************************** $0 $1 fail **********************"
    echo
fi
