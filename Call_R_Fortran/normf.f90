subroutine test_random(x, y) 
real*8 normrnd, unifrnd, x, y
call rndstart() 
x = normrnd()
y = unifrnd() 
call rndend() 
return
end