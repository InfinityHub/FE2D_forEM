MODULE FORTRAN_string_manipulation
  ! Procedures that manipulate strings.



  IMPLICIT NONE

  PRIVATE

  PUBLIC :: String_split_into_code_and_value, string_convert_to_lowercase

CONTAINS

  SUBROUTINE String_split_into_code_and_value( nLen, sLine, sCode, sValue, isComment )

    ! Parse a line read from a file into a code (name) & value.
    ! Force the code to be all lowercase with no ending colon.  Terminate
    ! the line at a '%' or '!' sign (these allow for user comments!)

    ! On return:
    ! Both code and value are strings which may contain spaces
    
    ! Args
    INTEGER, INTENT(IN)   :: nLen
    CHARACTER(nLen)       :: sLine
    CHARACTER(nlen), INTENT(OUT) :: sCode, sValue
    LOGICAL, INTENT(OUT)    :: isComment  ! whether the line is a comment line

    ! Local vars
    INTEGER :: iFrom, iTo

    ! Init returns
    isComment = .FALSE.
    sCode = ' '
    sValue = ' '

    ! Convert all tab characters to spaces (ichar returns the 'internal code number' of a character)
    FORALL( iTo = 1:nLen, ICHAR(sLine(iTo:iTo)) == 9 ) sLine(iTo:iTo) = ' '

    ! Skip any beginning blanks, iFrom presumably holds the position of the first non-blank char.
    DO iFrom = 1,nLen
       IF (sLine(iFrom:iFrom) /= ' ') EXIT
    ENDDO
    ! If the first char is a comment char, then the whole line is a comment.
    ! DGM April 2008 Also, if the line is blank, consider it a comment.
    IF (iFrom >= nLen) THEN !KWK may 2009 pulled this out in from since sometimes iFrom > nLen and this kills (iFrom:iFrom) below
       isComment = .TRUE.
       RETURN
    ENDIF

    ! A non-blank comment is defined as any lines starting with %, ! or #
    IF(  sLine(iFrom:iFrom) == '%' &
         .OR. sLine(iFrom:iFrom) == '!' &
         .OR. sLine(iFrom:iFrom) == '#') THEN
       isComment = .TRUE.
       RETURN
    ENDIF

    ! Pull off the code value. Cvt to lowercase as we go. (INDEX returns 0 if the given inquiry char is not found)
    iTo = INDEX(sLine,':') - 1
    IF (iTo < iFrom) THEN
       WRITE(*,*) 'Parsing Error: missing colon in line below:'
       WRITE(*,*) sLine
       RETURN
    ENDIF
    
    sCode = sLine(iFrom:iTo)
    CALL String_convert_to_Lowercase(sCode)

    ! Now, get the value after code
    ! Skip spaces after the colon; iFrom presumably holds the index position of the first non-blank char after semicolon
    DO iFrom = iTo+2, nLen
       IF (sLine(iFrom:iFrom) /= ' ') EXIT
    ENDDO

    ! Get the rest, up to any comment
    sValue = sLine(iFrom : nLen)
    iTo = LEN_TRIM(sValue)

    ! replace everything after '%' or '!' with spaces(blanks)
    iFrom = INDEX(sValue,'%')
    IF (iFrom > 0 .AND. iFrom < iTo) THEN
       sValue(iFrom:iTo) = ' '
    ENDIF
    
    iFrom = INDEX(sValue,'!')
    IF (iFrom > 0 .AND. iFrom < iTo) THEN
       sValue(iFrom:iTo) = ' '
    ENDIF
    !call Lower(sValue)   ! No: Some values are filenames which are case-sensitive on UNIX!

    RETURN
  end SUBROUTINE String_split_into_code_and_value

!==============================================================================!
  SUBROUTINE String_convert_to_Lowercase( s )
    IMPLICIT NONE
    CHARACTER(*), INTENT(INOUT)  :: s
    INTEGER i

    ! len_trim: returns the length of a character string ignoring any trailing blanks.
    DO  i=1, LEN_TRIM(s)
       IF  ( s(i:i) >= 'A' .AND. s(i:i) <= 'Z' ) THEN
          s(i:i) = CHAR(ICHAR(s(i:i)) + 32)
       ENDIF
    ENDDO

  END SUBROUTINE String_convert_to_Lowercase

end MODULE FORTRAN_string_manipulation
