program oler_hween;

//uses wincrt;

label start;

const maxn=10; maxk=3;


var a,b,c,d,e,f,g,h,u,v:real;
    i,j,k,n:integer;
    x:array[0..maxn] of real;
    y:array[0..maxn,0..maxk] of real;
    L:0..3;


function fxy (x,y:real):real;
begin
fxy:=a*x*sqr(x)+b*y*sqr(y)+c*sqr(x)+d*sqr(y)+e*x+f*y+g;
end;




Begin
start:writeln('        solving differential equation using Oler-Hween method.');
writeln('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
writeln('y`=f(x,y)');
writeln('what kind of equations do you want to solve?');
writeln('1:linear');writeln('2:second degree');writeln('3:third degree');writeln('0:get out');
readln(L);
a:=0;b:=0;c:=0;d:=0;e:=0;f:=0;g:=0;
 case l of
    1:begin
      writeln('y`=aX+bY+c');
      writeln('Input a,b,c');
      write('a=');read(e);write('b=');read(f);write('c=');read(g);
     end;
   2:begin
          writeln('y`=aX^2+bY^2+cX+dY+e');
          writeln('Input a,b,c,d,e please');
          write('a=');read(c);write('b=');read(d);write('c=');read(e);write('d=');read(f);write('e=');read(g);


     end;
   3:begin
          writeln('y`=aX^3+bY^3+cX^2+dY^2+eX+fY+g');
          writeln('Input a,b,c,d,e,f,g please');
          write('a=');read(a);write('b=');read(b);write('c=');read(c);
          write('d=');read(d);write('e=');read(e);write('f=');read(f);write('g=');read(g);

     end;
  0:halt;
 end;
writeln('X is in [u,v].Input u,v please');
write('u=');read(u);write('v=');read(v);
writeln('input number of steps n');
readln(n);
h:=(v-u)/n;
for i:=0 to n do
 x[i]:=u+i*h;
writeln('Input Y(X0) please');
readln(y[0,0]);
writeln('Finally please input number of improving steps k');
readln(k);
y[0,k]:=y[0,0];
for i:=0 to n-1 do
begin
 y[i+1,0]:=y[i,k]+h*fxy(x[i],y[i,k]);
 for j:=0 to k do
  y[i+1,j+1]:=y[i,k]+(h/2)*(fxy(x[i],y[i,k])+fxy(x[i+1],y[i+1,j]));

 writeln('y',i+1,'=',y[i+1,k+1]:0:9);
end;
goto start;
End.