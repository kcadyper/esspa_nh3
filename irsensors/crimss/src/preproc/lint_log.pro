function  lint_log, x1,x2, p1,p2, p0 
    return, x1*(p0/p1)^( alog(x2/x1)/alog(p2/p1) )
end
