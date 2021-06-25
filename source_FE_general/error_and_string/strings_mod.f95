! file I/O
MODULE strings_mod
! Contains various subroutines for working with strings.
   
   USE error_cls
   
   IMPLICIT NONE
   
   PRIVATE
   PUBLIC :: PATHLEN, WORDLEN, PRINTLEN
   PUBLIC :: compare_strings, remove_comments, replace_char, replace_substring, replace_commas, replace_spaces, all_spaces, &
             file_parts, rel2abs, concatenate, time_string, tokenize, count_words

!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
! STATIC PARAMETERS
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
   
   INTEGER, PARAMETER :: PATHLEN=256 ! length of character strings that hold full path names to files
   INTEGER, PARAMETER :: WORDLEN=PATHLEN ! length of character strings that hold a tokenized word
   INTEGER, PARAMETER :: PRINTLEN=128 ! length of character strings that hold something to print to screen

CONTAINS
   
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
! COMPARISON
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

   LOGICAL PURE FUNCTION compare_strings(s1,s2) RESULT(com)
   ! Returns true if two strings are the same.
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: s1,s2
      INTEGER :: n1,n2
      n1 = LEN_TRIM(s1)
      n2 = LEN_TRIM(s2)
      IF (n1==n2) THEN
         com = ( s1(1:n1) == s2(1:n2) )
      ELSE
         com = .FALSE.
      END IF
   END FUNCTION compare_strings

   
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
! REMOVAL AND REPLACEMENT
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

   PURE SUBROUTINE remove_comments(tt,cc)
   ! Removes comments (anything at or after first comment character) from a character string.
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(INOUT) :: tt ! character string to remove comments from
      CHARACTER(LEN=1), INTENT(IN) :: cc ! comment character
      INTEGER :: ind,i
      ! Find first comment character:
      ind = INDEX(tt,cc,.false.) ! Returns the first starting position for a substring within a string.
      IF (ind<=0) RETURN ! comment character not found
      ! Remove comment:
      DO i=ind,LEN(tt)
         tt(i:i) = ' ' ! replace with spaces
      END DO
   END SUBROUTINE remove_comments

!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

   PURE SUBROUTINE replace_char(t,c1,c2)
   ! Replaces any occurrences of the character c1 in string t with character c2.
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(INOUT) :: t ! character string to use
      CHARACTER(LEN=1), INTENT(IN) :: c1,c2 ! characters to look for and replace by
      INTEGER :: i
      ! Loop over each character in the string:
      DO i=1,LEN(t)
         ! Check for c1 and replace with c2 if found:
         IF (t(i:i)==c1(1:1)) t(i:i)=c2
      END DO
   END SUBROUTINE replace_char

!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

   PURE SUBROUTINE replace_substring(t,s1,s2,error)
   ! Replaces any occurrences of the substring s1 in string t with substring c2.
   ! The two substrings must be equal length.
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(INOUT) :: t ! character string to use
      CHARACTER(LEN=*), INTENT(IN) :: s1,s2 ! substrings to look for and replace by
      TYPE(error_type), INTENT(INOUT) :: error
      INTEGER :: n,i
      ! Check inputs:
      n = LEN(s1)
      IF (LEN(s2)/=n) THEN
         CALL error_construct(error,ERROR_GENERAL,'strings_mod','replace_substring','substrings inconsistent lengths')
         RETURN
      END IF
      IF (compare_strings(s1,s2)) THEN
         CALL error_construct(error,ERROR_GENERAL,'strings_mod','replace_substring','substrings are equal')
         RETURN
      END IF
      n = n - 1
      ! Loop until all substrings removed:
      i = INDEX(t,s1) ! find first substring s1
      DO WHILE(i>0) ! while substring s1 exists in t
         t(i:(i+n)) = s2 ! replace substring s1 with s2
         i = INDEX(t,s1) ! find next substring s1
      END DO
   END SUBROUTINE replace_substring

!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

   PURE SUBROUTINE replace_commas(t)
   ! Replaces any commas in string t with spaces.
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(INOUT) :: t ! character string to use
      CALL replace_char(t,',',' ')
   END SUBROUTINE replace_commas

!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

   PURE SUBROUTINE replace_spaces(t,c)
   ! Replaces spaces in character string t with the character c.
   ! Trailing blanks are not altered.
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(INOUT) :: t ! character string to replace spaces in
      CHARACTER(LEN=1), INTENT(IN) :: c ! character to replace spaces by
      INTEGER :: i
      ! Loop over each character in the string except for trailing blanks:
      DO i=1,LEN_TRIM(t)
         ! Check for a space and replace if found:
         IF (t(i:i)==' ') t(i:i)=c
      END DO
   END SUBROUTINE replace_spaces
   
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

   PURE SUBROUTINE all_spaces(t)
   ! Sets all characters in string t with spaces.
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(INOUT) :: t
      INTEGER :: i
      DO i=1,LEN_TRIM(t)
         t(i:i) = ' '
      END DO
   END SUBROUTINE all_spaces
   
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
! STRINGS THAT REPRESENT FILE NAMES
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

   PURE SUBROUTINE file_parts(filename,path,name,ext)
   ! Extracts file path, name and extension from an full input file name.
   ! The extension includes the dot character.
   ! The path includes the final / or \ character.
   
   ! There is some test code for this subroutine at the top of vinv.f95.
      
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: filename
      CHARACTER(LEN=*), INTENT(OUT) :: path, name
      CHARACTER(LEN=*), INTENT(OUT) :: ext
      
      INTEGER :: iext,iname
      INTEGER :: flen ! length of filename
      
      ! Initialize output strings to all spaces (just in case):
      path = ' '
      name = ' '
      ext = ' '
      
      ! Find file part information:
      CALL find_parts(filename,iname,iext,flen)
      
      ! Set output strings:
      IF (iname>1) path = filename(1:(iname-1)) ! includes the last '/' or '\' character
      IF (iext<=0) THEN ! '.' not found
         IF (iname>=1) name = filename(iname:flen)
         ext = ' '
      ELSE
         IF (iname>=1) name = filename(iname:(iext-1))
         ext = filename(iext:flen) ! includes the dot character
      END IF
      
   END SUBROUTINE file_parts

!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

   PURE SUBROUTINE find_parts(filename,iname,iext,flen)
   ! Finds indices marking the start of different parts of a file name.
   !
   ! iname is the index after the last '/' or '\' character
   ! iname=1 if '/' or '\' not found
   ! iname=0 under some special circumstances (e.g. filename='..')
   !
   ! iext is the index of the '.' character
   ! iext=0 if '.' not found
   !
   ! flen = length of filename without any trailing blanks (result of LEN_TRIM)
   
   ! There is some test code for this subroutine at the top of vinv.f95.
   
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: filename
      INTEGER, INTENT(OUT) :: iname, iext
      INTEGER, INTENT(OUT)  :: flen ! length of filename
      
      INTEGER :: i1,i2
      
      ! Determine length of input file name:
      flen = LEN_TRIM(filename) ! Returns the length of the string without the possibly trailing blanks.
      
      ! Check for the only case I've found that doesn't work with the logic below:
      IF (filename(1:flen)=='..') THEN
         iname = 0
         iext = 0
      END IF
      
      ! Find the last occurrence of '\' or '/':
      i1 = INDEX(filename,'/',.TRUE.)
      i2 = INDEX(filename,'\',.TRUE.)
      iname = MAX(i1,i2)
      iname = iname + 1
      iname = MAX(iname,1) ! =1 if '\' or '\' not found
      
      ! Find extension:
      iext = INDEX(filename(iname:flen),'.',.TRUE.) ! Returns the last starting position for a substring within a string.
      IF (iext<=0) THEN
         iext = 0 ! =0 if '.' not found
      ELSE
         iext = iname + iext - 1
      END IF
      
   END SUBROUTINE find_parts

!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

   SUBROUTINE rel2abs(path1,path2,error)
   ! Converts path2 to absolute path based on path1.
   ! The path1 input should be any path (possibly including file name at end).
   ! Path2 can be an absolute path or some path relative to path1.
   ! For example, with input parameters
   !   path1 = /Users/Peter/Documents/temp/name.ext
   !   path2 = ../../Work/inv/doc.txt
   ! the outputs will be
   !   path1 = /Users/Peter/Documents/temp/name.ext (unchanged)
   !   path2 = /Users/Peter/Work/inv/doc.txt (changed)
   !
   ! Be careful: if path1 doesn't end in / (or \) then it may be truncated (refer to subroutine file_parts).
   !
   ! If path1 is empty then path2 will not change.
   ! If path2 is an absolute path (staring with /,\,~) then it is not altered.
   ! If path2 is a forward relative path (not starting with .,/,\,~) or simply a file name then path1 is prepended onto path2.
   ! If path2 is empty or = 'null' then path2 will not change.
   !
   ! The path1 input is never altered (has intent(in)).
   
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: path1
      CHARACTER(LEN=*), INTENT(INOUT) :: path2
      TYPE(error_type), INTENT(INOUT) :: error
      
      INTEGER :: i1,i2,i,plen
      CHARACTER(LEN=PATHLEN) :: p1,p2,p,n,e
      
      ! Get file parts for path1:
      CALL file_parts(path1,p1,n,e) ! strips off any name.ext; p1 ends in / or \
      
      ! Make sure path2 doesn't start with blanks:
      p2 = ADJUSTL(path2)
      
      ! Make sure path1 doesn't start with blanks:
      IF (path1(1:1)==' ') THEN
            CALL error_construct(error,ERROR_GENERAL,'strings_mod','rel2abs','path1 must not start with blanks')
            RETURN
      END IF
      ! Make sure path1 isn't a relative path:
      IF (path1(1:1)=='.') THEN
            CALL error_construct(error,ERROR_GENERAL,'strings_mod','rel2abs','path1 must not be a relative path')
            RETURN
      END IF
      
      ! Check for special cases:
      IF (LEN_TRIM(p1)==0) RETURN ! path1 is empty
      IF (LEN_TRIM(p2)==0) RETURN ! path2 is empty
      IF (p2(1:4)=='null') RETURN ! special 'null' value
      IF ((p2(1:1)=='/').OR.(p2(1:1)=='\').OR.(p2(1:1)=='~')) RETURN ! assume path2 is absolute path
      IF (p2(1:1)/='.') THEN ! path2 might be a forward relative path or simply a file name
         CALL concatenate(error,path2,p1,p2)
         IF (error_check(error)) CALL error_pass(error,'strings_mod','rel2abs')
         RETURN
      END IF
      
      ! Loop until finished:
      DO WHILE(.TRUE.)
         
         ! Check for relative path name:
         IF (p2(1:1)/='.') EXIT ! exits once all '../' have been stripped off
         
         ! Find the first occurrence of '..\' or '../' in the relative path:
         i1 = INDEX(p2,'../')
         i2 = INDEX(p2,'..\')
         IF ((i1>0).AND.(i2>0)) THEN
            CALL error_construct(error,ERROR_GENERAL,'strings_mod','rel2abs','mixed ..\ and ../')
            RETURN
         END IF
         i = MAX(i1,i2) ! ONE is <=0
         IF (i<=0) THEN
            CALL error_construct(error,ERROR_GENERAL,'strings_mod','rel2abs','something wrong: perhaps .\ or ./')
            RETURN
         END IF
         
         ! Strip off the first '..\' or '../' from the relative path:
         i = i + 3
         plen = LEN_TRIM(p2)
         p2 = p2(i:plen)
         
         ! Find the last occurrence of '\' or '/' in the first input path:
         i1 = INDEX(p1,'/',.TRUE.)
         i2 = INDEX(p1,'\',.TRUE.)
         i = MAX(i1,i2)
         IF (i<=0) THEN
            CALL error_construct(error,ERROR_GENERAL,'strings_mod','rel2abs','unbalanced paths (i)')
            RETURN
         END IF
         
         ! Strip off the last part of the path in the first input path:
         plen = LEN_TRIM(p1) - 1 ! -1 to remove the '/' or '\' character
         p = p1(1:plen)
         CALL file_parts(p,p1,n,e)
         plen = LEN_TRIM(p1)
         IF (plen<=0) THEN
            PRINT *, plen
            PRINT *, 'path1: "',TRIM(path1),'"'
            PRINT *, 'path2: "',TRIM(path2),'"'
            PRINT *, 'p1   : "',TRIM(p1),'"'
            PRINT *, 'p2   : "',TRIM(p2),'"'
            CALL error_construct(error,ERROR_GENERAL,'strings_mod','rel2abs','unbalanced paths (plen)')
            RETURN
         END IF
         
      END DO

      ! Combine the two paths:
!      PRINT *, 'path1: "',TRIM(path1),'"'
!      PRINT *, 'path2: "',TRIM(path2),'"'
!      PRINT *, 'p1   : "',TRIM(p1),'"'
!      PRINT *, 'p2   : "',TRIM(p2),'"'
      CALL concatenate(error,path2,p1,p2)
      IF (error_check(error)) CALL error_pass(error,'strings_mod','rel2abs')
!      PRINT *, 'path2: "',TRIM(path2),'"'
!      PRINT *, ' '
      
   END SUBROUTINE rel2abs

!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
! COMBINING, CONSTRUCTING
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

   PURE SUBROUTINE concatenate(error,s,s1,s2,s3_op,s4_op,trimflag_op)
   ! Concatenates two or three strings.
   ! Checks that the output string is large enough to hold the trimmed contents of the other strings.
   ! If any of the strings have trailing spaces then those spaces will get removed unless trimflag=.FALSE.
   ! Leading spaces are never removed.
   ! Avoids use of the // operator.
      
      IMPLICIT NONE
      TYPE(error_type), INTENT(INOUT) :: error
      CHARACTER(LEN=*), INTENT(OUT) :: s
      CHARACTER(LEN=*), INTENT(IN) :: s1,s2
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: s3_op,s4_op
      LOGICAL, INTENT(IN), OPTIONAL :: trimflag_op
      
      LOGICAL :: present_s3, present_s4, trimflag
      INTEGER :: n1,n2,n3,n4,n,i
      
      ! Check optional inputs:
      present_s3 = PRESENT(s3_op)
      present_s4 = PRESENT(s4_op)
      IF (PRESENT(trimflag_op)) THEN
         trimflag = trimflag_op
      ELSE
         trimflag = .TRUE.
      END IF
      
      ! Get trimmed length of the strings:
      IF (trimflag) THEN
         n1 = LEN_TRIM(s1)
         n2 = LEN_TRIM(s2)
         IF (present_s3) THEN
            n3 = LEN_TRIM(s3_op)
            IF (present_s4) THEN
               n4 = LEN_TRIM(s4_op)
            ELSE
               n4 = 0
            END IF
         ELSE
            n3 = 0
            n4 = 0
         END IF
      ELSE
         n1 = LEN(s1)
         IF (present_s3) THEN
            n2 = LEN(s2)
            IF (present_s4) THEN
               n3 = LEN(s3_op)
               n4 = LEN_TRIM(s4_op) ! n4 is the last string in this case
            ELSE
               n3 = LEN_TRIM(s3_op) ! n3 is the last string in this case
               n4 = 0
            END IF
         ELSE
            n2 = LEN_TRIM(s2) ! n2 is the last string in this case
            n3 = 0
            n4 = 0
         END IF
      END IF
      n = LEN(s) ! amount of space available for concatenation
      
      ! Check output string is long enough:
      IF ( n < (n1+n2+n3+n4) ) THEN
         CALL error_debug(error,'strings_mod','concatenate','output string too short:',(/n,(n1+n2+n3+n4)/))
         RETURN
      END IF
      
      ! Clear the string:
      DO i=1,n
         s(i:i) = ' '
      END DO
      
      ! Concatenate:
      s(1:n1) = s1(1:n1)
      s((n1+1):(n1+n2)) = s2(1:n2)
      IF (present_s3) THEN
         s((n1+n2+1):(n1+n2+n3)) = s3_op(1:n3)
         IF (present_s4) s((n1+n2+n3+1):(n1+n2+n3+n4)) = s4_op(1:n4)
      END IF
      
   END SUBROUTINE concatenate
   
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

   SUBROUTINE time_string(t,s,error)
   ! Puts the wall time in t into the string s
   ! The format is
   !    year/month/day ; hour:minute:second
   ! with this many characters
   !    yyyy/mm/dd ; hh:mm:ss
      IMPLICIT NONE
      INTEGER, DIMENSION(8), INTENT(IN) :: t ! the values output from DATE_AND_TIME
      CHARACTER(LEN=*), INTENT(OUT) :: s
      TYPE(error_type), INTENT(INOUT) :: error
      IF (LEN(s)<21) THEN
         CALL error_construct(error,ERROR_GENERAL,'numerics_mod','time_string','not enough space in string (need len>=21)')
         RETURN
      END IF
      WRITE(UNIT=s,FMT='(I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)') &
         t(1),'/',t(2),'/',t(3),' ; ',t(5),':',t(6),':',t(7)
   END SUBROUTINE time_string

!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
! SPLITTING
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

   SUBROUTINE tokenize(str,tok,tmtao,nmax,n,words,error)
   ! Tokenizes a string.
   ! If tok is a space then any tab characters are treated as spaces.
      
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: str ! string to tokenize
      CHARACTER(LEN=1), INTENT(IN) :: tok ! token character
      LOGICAL, INTENT(IN) :: tmtao ! treat multiple tokens as ONE?
      INTEGER, INTENT(IN) :: nmax ! number of words to read (set to <=0 to read all possible words)
      INTEGER, INTENT(OUT) :: n ! number of words found
      CHARACTER(LEN=*), DIMENSION(:), INTENT(OUT) :: words ! output words
      TYPE(error_type), INTENT(INOUT) :: error
      
      INTEGER :: i,pos,leng
      CHARACTER(LEN=WORDLEN) :: wd
      CHARACTER(LEN=LEN(str)) :: s
      LOGICAL :: ok
      
      ! Get length of string:
      leng = LEN_TRIM(str)
      
      ! If tok is a space then replace all tab characters with spaces:
      s = ' '
      s(1:leng) = str(1:leng)
      IF (tok(1:1)==' ') THEN
         DO i=1,leng
            IF (IACHAR(s(i:i))==9) THEN
               s(i:i) = ' '
            END IF
         END DO
      END IF
      
      ! Check for trivial solution:
      n = 0 ! counts the number of words read
      IF (leng<1) RETURN
      
      ! Initialize all words to blanks:
      DO i=1,SIZE(words,1)
         words(i) = ' '
      END DO
      
      ! Loop until reached end of the string:
      pos = 1
      DO
         ! Get the next word:
         CALL next_word(s,tok,tmtao,pos,wd,ok,error)
         IF (error_check(error)) THEN
            CALL error_pass(error,'strings_mod','tokenize')
            RETURN
         END IF
         ! Check for reaching the end of the string
         IF (.NOT.ok) EXIT
         ! Increment the word counter:
         n = n + 1
         ! Check output array is long enough:
         IF (n>SIZE(words,1)) THEN
            PRINT *, n, SIZE(words,1)
            PRINT *, 'str = "',TRIM(str),'"'
            PRINT *, 's   = "',TRIM(s),'"'
            PRINT *, 'tok = "',tok,'"'
            PRINT *, tmtao,nmax
            CALL error_construct(error,ERROR_GENERAL,'strings_mod','tokenize','not enough room in words array')
            RETURN
         END IF
         ! Insert word into output array:
         words(n) = wd
         ! Check for reading enough words:
         IF ((nmax>0).AND.(n>=nmax)) EXIT
      END DO
      
      ! Check if we read the correct number of words:
      IF ((nmax>0).AND.(n/=nmax)) THEN
         PRINT *, n, SIZE(words,1)
         PRINT *, 'str = "',TRIM(str),'"'
         PRINT *, 's   = "',TRIM(s),'"'
         PRINT *, 'tok = "',tok,'"'
         PRINT *, 'IACHAR(tok) = ',IACHAR(tok)
         PRINT *, tmtao,nmax
         CALL error_construct(error,ERROR_GENERAL,'strings_mod','tokenize','not enough words found')
         RETURN
      END IF
      
   END SUBROUTINE tokenize
   
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

   SUBROUTINE count_words(str,tok,tmtao,n,error)
   ! Counts how many words are in a string.
      
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: str ! string to tokenize
      CHARACTER(LEN=1), INTENT(IN) :: tok ! token character
      LOGICAL, INTENT(IN) :: tmtao ! treat multiple tokens as ONE?
      INTEGER, INTENT(OUT) :: n ! the number of words
      TYPE(error_type), INTENT(INOUT) :: error
      
      INTEGER :: pos,leng
      CHARACTER(LEN=WORDLEN) :: wd
      LOGICAL :: ok
      
      ! Get length of string:
      leng = LEN_TRIM(str)
      !PRINT *, '"',TRIM(str),'"'
      
      ! Check for trivial solution:
      n = 0 ! counts the number of words read
      IF (leng<1) RETURN
      
      ! Loop until reached end of the string:
      pos = 1
      DO
         ! Get the next word:
         CALL next_word(str,tok,tmtao,pos,wd,ok,error)
         IF (error_check(error)) THEN
            CALL error_pass(error,'strings_mod','count_words')
            RETURN
         END IF
         ! Check for reaching the end of the string
         IF (.NOT.ok) EXIT
         ! Increment the word counter:
         n = n + 1
         !PRINT *, n,pos,'"',TRIM(wd),'"'
      END DO
      
   END SUBROUTINE count_words
   
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

   SUBROUTINE next_word(str,tok,tmtao,pos1,word,ok,error)
   ! Gets next word in a string with tokens.
   ! Returns true if a word was found, false if reached the end of the string.
      
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: str ! string to tokenize
      CHARACTER(LEN=1), INTENT(IN) :: tok ! token character
      LOGICAL, INTENT(IN) :: tmtao ! treat multiple tokens as ONE?
      INTEGER, INTENT(INOUT) :: pos1 ! position to start from; on exit is the position to start looking for the next word
      CHARACTER(LEN=*), INTENT(OUT) :: word ! output word
      LOGICAL, INTENT(OUT) :: ok
      TYPE(error_type), INTENT(INOUT) :: error
      
      INTEGER :: pos2,strleng,wordleng
      
      ! Initialize ok flag:
      ok = .FALSE.
      
      ! Get length of string:
      strleng = LEN(str)
      
      ! Check for pos1 out of range:
      IF (pos1<=0) THEN
         CALL error_construct(error,ERROR_GENERAL,'strings_mod','next_word','pos1 too low')
         RETURN
      END IF
      
      ! Find end of word:
      DO ! this loop is required in case tmtao=.TRUE.
         
         ! Check for reaching the end of the string:
         IF (pos1>strleng) RETURN ! ok=.FALSE.
         
         ! Find next index of the token character:
         pos2 = INDEX(str(pos1:strleng),tok)
         
         ! Adjust value of pos2 as required:
         IF (pos2==0) THEN ! no token character found beyond the starting position
            pos2 = strleng ! the next word extends to the end of the string
         ELSE IF (pos2>0) THEN ! token character found
            pos2 = pos1 + pos2 - 2 ! the index before the next token character
         ELSE
            CALL error_construct(error,ERROR_GENERAL,'strings_mod','next_word','unexpected value of pos2')
            RETURN
         END IF
         wordleng = pos2 - pos1 + 1 ! length of the word
         
         ! Watch for the case where there are multiple tokens in a row:
         IF (tmtao.AND.(wordleng==0)) THEN
            ! Reset pos1 and cycle:
            pos1 = pos1 + 1
            CYCLE
         END IF
         
         ! Extract the word, watching for an empty word:
         IF (wordleng>0) THEN ! a word of ONE or more characters was found
            word = str(pos1:pos2)
         ELSE IF (wordleng==0) THEN ! the word is empty
            word = ' '
         ELSE
            CALL error_construct(error,ERROR_GENERAL,'strings_mod','next_word','unexpected value of pos2-pos1')
            RETURN
         END IF
         
         ! Alter the pos1 value to point to the position in the string that starts the next word:
         pos1 = pos2 + 2
         
         EXIT
      END DO
      
      ! Set ok flag to true:
      ok = .TRUE.
      
   END SUBROUTINE next_word
   
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

END MODULE strings_mod
