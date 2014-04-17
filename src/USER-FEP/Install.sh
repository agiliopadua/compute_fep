# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

# this is default Install.sh for all packages
# if package has an auxiliary library or a file with a dependency,
# then package dir has its own customized Install.sh

mode=$1

# arg1 = file, arg2 = file it depends on

action () {
  if (test $mode = 0) then
    rm -f ../$1
  elif (! cmp -s $1 ../$1) then
    if (test -z "$2" || test -e ../$2) then
      cp $1 ..
      if (test $mode = 2) then
        echo "  updating src/$1"
      fi
    fi
  elif (test -n "$2") then
    if (test ! -e ../$2) then
      rm -f ../$1
    fi
  fi
}

## backup fix_adapt
#if (test $mode != 0) then
#  if (test -e ../fix_adapt.cpp) then
#    if (test ! -d ../ORIG.fix_adapt) then
#      mkdir ../ORIG.fix_adapt
#      mv ../fix_adapt.* ../ORIG.fix_adapt
#      echo "  backed up fix_adapt sources to ORIG.fix_adapt directory"
#    else
#      echo "  backup orig.fix_adapt directory present, not overwritten"
#    fi
#  fi
#fi

# all package files with no dependencies

for file in *.cpp *.h; do
  action $file
done

## restore fix adapt
#if (test $mode = 0) then
#  if (test -d ../ORIG.fix_adapt) then
#    mv -f ../ORIG.fix_adapt/* ..
#    rmdir ../ORIG.fix_adapt
#    echo "  restored fix_adapt"
#  fi
#fi
