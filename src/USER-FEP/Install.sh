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

# backup fix_adapt
if (test $mode != 0) then
  if (test -e ../fix_adapt.cpp) then
    if (test ! -e orig.fix_adapt.cpp) then
      mv ../fix_adapt.h orig.fix_adapt.h
      mv ../fix_adapt.cpp orig.fix_adapt.cpp
      echo "  backed up fix_adapt to USER-FEP/orig.fix_adapt"
    else
      echo "  backup USER-FEP/orig.fix_adapt present, not overwritten"
    fi
  fi
fi

# all package files with no dependencies

for file in *.cpp *.h; do
  action $file
done

# restore fix adapt
if (test $mode = 0) then
  if (test -e orig.fix_adapt.cpp) then
    mv orig.fix_adapt.h ../fix_adapt.h
    mv orig.fix_adapt.cpp ../fix_adapt.cpp
    echo "  restored fix_adapt"
  fi
fi
