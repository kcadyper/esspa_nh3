pro run_gfs, other_args

compile_opt strictarr

  args = command_line_args(count=nargs)
  if nargs gt 0L then begin
    file=args[0]
    outpath='./'
    if nargs gt 1L then outpath=args[1] 
    convert_gfs_grib2, file, outpath=outpath
  endif
    
  if keyword_set(other_args) then begin
    print, 'outpath:', outpath
    help, other_args
    if (n_elements(other_args) gt 0L) then print, other_args
  endif  
end
