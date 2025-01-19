program newton_gaus;

uses wincrt;

label start;
const maxn=3;

type matrix=array[1..maxn,1..maxn] of real;
     vector=array[1..maxn] of real;

var  Jak:matrix;
     U,S,DS,prev:vector;
     a,b,c,d,e,f,g,h,i,j:array[1..3] of integer;
     k:1..3;
     error:real;
     st,L:byte;
function Fxyz(x,y,z:real;k:byte):real;
begin
 Fxyz:=a[k]*x*sqr(x)+b[k]*y*sqr(y)+c[k]*z*sqr(z)+d[k]*sqr(x)+e[k]*sqr(y)+f[k]*sqr(z)+g[k]*x+h[k]*y+i[k]*z+j[k];
end;

function derive(a,b,c:integer;x:real):real;
begin
 derive:=3*a*sqr(x)+2*b*x+c;
end;

Procedure Gaus(oparray:matrix;frvec:vector;n:byte;var X:vector);

type expanded=array[1..maxn,1..maxn+1,0..maxn] of real;
     GausArray=array[1..maxn,1..maxn+1] of real;

var  A:expanded; B:GausArray;
     i,j,k:0..maxn;
     tempsum:real;

begin

 for i:=1 to n do
 begin

  for j:=1 to n do
   A[i,j,0]:=OpArray[i,j];
  A[i,n+1,0]:=FrVec[i];

 end;

 for k:=1 to n do
 begin

  for j:=k to n+1 do
   B[k,j]:=A[k,j,k-1]/A[k,k,k-1];

  for i:=k+1 to n do
   for j:=k to n+1 do
      a[i,j,k]:=a[i,j,k-1]-a[i,k,k-1]*B[k,j];

 end;

 X[n]:=B[n,n+1];
 for i:=n-1 downto 1 do
 begin
  tempsum:=0;
  for j:=i+1 to n do
   tempsum:=tempsum+X[j]*B[i,j];
  X[i]:=B[i,n+1]-tempsum;
 end;
end;

BEGIN {main program}
 start:writeln('        solving nonlinear equations using advanced Newton-Gauss method.');
 writeln('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
 for k:=1 to maxn do
 begin
 a[k]:=0;b[k]:=0;c[k]:=0;d[k]:=0;e[k]:=0;f[k]:=0;g[k]:=0;h[k]:=0;i[k]:=0;j[k]:=0;
 end;
 writeln('what kind of equations do you want to solve?');
 writeln('1:linear');writeln('2:second degree');writeln('3:third degree');writeln('0:get out');
 readln(L);
 case l of
   1:begin
      for k:=1 to maxn do
       begin
         writeln('input the operators of U',k,'=aX+bY+cZ+d');
         write('a',k,'=');read(jak[k,1]);write('b',k,'=');read(jak[k,2]);
         write('c',k,'=');read(jak[k,3]);write('d',k,'=');read(u[k]);
       end;
      Gaus(Jak,U,3,DS);
      writeln('X=',DS[1]:0:7);
      writeln('Y=',DS[2]:0:7);
      writeln('Z=',DS[3]:0:7);
      goto start;
     end;
     
   2:for k:=1 to maxn do
     begin
         writeln('input the operators of U',k,'=aX^2+bY^2+cZ^2+dX+eY+fZ+g');
         write('a',k,'=');read(d[k]);write('b',k,'=');read(e[k]);write('c',k,'=');read(f[k]);write('d',k,'=');
         read(g[k]);write('e',k,'=');read(h[k]);write('f',k,'=');read(i[k]);write('g',k,'=');read(j[k]);
     end;


   3:for k:=1 to maxn do
     begin
        writeln('input the operators of U',k,'=aX^3+bY^3+cZ^3+dX^2+eY^2+fZ^2+gX+hY+iZ+j');
        write('a',k,'=');read(a[k]);write('b',k,'=');read(b[k]);write('c',k,'=');read(c[k]);write('d',k,'=');read(d[k]);
        write('e',k,'=');read(e[k]);write('f',k,'=');read(f[k]);write('g',k,'=');read(g[k]);write('h',k,'=');read(h[k]);
        write('i',k,'=');read(i[k]);write('j',k,'=');read(j[k]);
     end;
   4:for k:=1 to maxn do
     begin
        writeln('input the operators of U',k,'=aX^3+bY^3+cZ^3+dX^2+eY^2+fZ^2+gX+hY+iZ+j');
        write('a',k,'=');read(a[k]);write('b',k,'=');read(b[k]);write('c',k,'=');read(c[k]);write('d',k,'=');read(d[k]);
        write('e',k,'=');read(e[k]);write('f',k,'=');read(f[k]);write('g',k,'=');read(g[k]);write('h',k,'=');read(h[k]);
        write('i',k,'=');read(i[k]);write('j',k,'=');read(j[k]);
     end;

   0:halt;
 end;
 readln;
 writeln(' please input first step X0,Y0,Z0 ');
 write('X',0,'=');read(s[1]);
 write('Y',0,'=');read(s[2]);
 write('Z',0,'=');read(s[3]);
 writeln('input error ratio please');
 readln(error);
 st:=1;
 repeat
  for k:=1 to 3 do
  begin
   U[k]:=-Fxyz(s[1],s[2],s[3],k);
   Jak[k,1]:=derive(a[k],d[k],g[k],s[1]);
   Jak[k,2]:=derive(b[k],e[k],h[k],s[2]);
   Jak[k,3]:=derive(c[k],f[k],i[k],s[3]);
  end;
  Gaus(Jak,U,3,DS);
  for k:=1 to 3 do
    begin
     prev[k]:=S[k];
     S[k]:=S[k]+DS[k];
    end;
   writeln('X',st,'=',S[1]:0:7);
   writeln('Y',st,'=',S[2]:0:7);
   writeln('Z',st,'=',S[3]:0:7);
   writeln;
   st:=st+1;
 until abs(s[1]-prev[1])/abs(s[1])<=error;
 goto start;
END.
