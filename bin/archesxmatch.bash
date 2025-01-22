#!/bin/bash
declare -rx host=serendib.unistra.fr #3:8080
declare -rx url=${host}/ARCHESWebService/XMatchARCHES
declare -rx cookie=.logcookie

usage(){
  echo ""
  echo "USAGE   : $0 ACTION [FILE]                                           "
  echo "ACTIONS :                                                            "
  echo " h, help            Print this usage message, then exit.             "
  echo " i, login           Log in the web service.                          "
  echo " q, quit            Log out the web service.                         "
  echo " c, chgpasswd       Change your password.                            "
  echo " ?, logged          Check if your are logged into the web service.   "
  echo " l, ls              List the files in  your distant repository.      "
  echo " p, put       FILE  Upload a file into your distant repository.      "
  echo " g, get       FILE  Download a file from your distant repository.    "
  echo " r, rm        FILE  Remove a file from your distant repository.      "
  echo " x, xmatch    FILE  Execute the given xmatch script file.            "
  echo ""
}

###########
## TOOLS ##

clean(){
  rm ${cookie}
}

checkfile(){
  [[ ! -f $1 ]] && { echo "'$1' is not a file!"; exit 1; }
}

# 1: password to be encoded
encodePasswd(){
  declare -rx a=$(echo -n "$1" | openssl dgst -sha1 -binary | xxd -p)
  echo -n "${a}" | openssl enc -a
}

#############
## ACTIONS ##

chgpasswd(){
  read -p "Enter you login: " login
  read -s -p "Enter old password: " oldPasswd; echo
  read -s -p "Enter new password: " newPasswd; echo
  read -s -p "Reenter new password: " newPasswdBis; echo
  # Verif both new passwd are the same
  [[ ${newPasswd} != ${newPasswdBis} ]] && { \
      echo "Operation aborted: both new password are different!"; \
      exit 0; \
  }
  # A form-based authentication to change password
  [[ -f ${cookie} ]] && { rm ${cookie}; }
  curl --junk-session-cookies \
       --data "cmd=chgpasswd&username=${login}&oldpassword=$(encodePasswd ${oldPasswd})&newpassword=$(encodePasswd ${newPasswd})" \
       ${url} 
}

login(){
  read -p "Enter you login: " login
  read -s -p "Enter Password: " passwd; echo
  # A form-based authentication
  curl --junk-session-cookies --cookie-jar ${cookie} \
       --data "cmd=login&username=${login}&password=$(encodePasswd ${passwd})" \
       ${url} 
}

islogged(){
  curl --cookie ${cookie} ${url}"?cmd=islogged"
}

quit(){
  curl --cookie ${cookie} ${url}"?cmd=quit"
  clean
}

list(){
  curl --cookie ${cookie} ${url}"?cmd=ls"
}

# 1: file to be uploaded
upload(){
  curl --cookie ${cookie} -X POST -F cmd=put -F "$1=@$1" ${url}
}

# 1: name of the file to be removed
remove(){
  curl --cookie ${cookie} --data "cmd=rm&fileName=$1" ${url}
}

# 1: xmatch script file
xmatch(){
  #curl --cookie ${cookie} -X POST -F cmd=xmatch -F "notimeout=true" -F "script=@$1" ${url}
  curl --cookie ${cookie} -X POST -F cmd=xmatch -F "script=@$1" ${url}
}

# 1: file to be downloaded
download(){
  curl --cookie ${cookie} --data "cmd=get&fileName=$1" -o "$1" ${url}
}

###########
## TESTS ##

tests(){
  $0 help
  $0 logged
  $0 login
  $0 logged
  $0 ls
  $0 get "file blabla"
  $0 rm "file blabla"
  $0 ls
  $0 put "file blabla"
  $0 ls
  $0 xmatch xmatchscript.txt
  $0 ls
  $0 quit
  $0 logged
  exit 0
}

#############
## PARSING ##

if [[ $# == 1 ]] ; then
  [[ ${USER} == 'pineau' && $1 == 'test' ]] && { tests; exit 0; }
  [[ $1 == '-h' || $1 == '--help' ]] && { usage; exit 0; }
  [[ $1 == 'h' || $1 == 'help' ]] && { usage; exit 0; }
  [[ $1 == 'i' || $1 == 'login' ]] && { login; exit 0; }
  [[ $1 == 'q' || $1 == 'quit' ]] && { quit; exit 0; }
  [[ $1 == 'c' || $1 == 'chgpasswd' ]] && { chgpasswd; exit 0; }
  [[ $1 == '?' || $1 == 'logged' ]] &&  { islogged; exit 0; }
  [[ $1 == 'l' || $1 == 'ls' ]] && { list; exit 0; }
elif [[ $# == 2 ]]; then
  [[ $1 == 'p' || $1 == 'put' ]] && { checkfile "$2"; upload "$2"; exit 0; }
  [[ $1 == 'g'  || $1 == 'get' ]] && { download "$2"; exit 0;}
  [[ $1 == 'r'  || $1 == 'rm' ]] && { remove "$2"; exit 0;}
  [[ $1 == 'x'  || $1 == 'xmatch' ]] && { checkfile "$2"; xmatch "$2"; exit 0;}
fi
usage; exit 1;

