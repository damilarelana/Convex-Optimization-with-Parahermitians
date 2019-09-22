
function [deltatoscalar,scalartodscales] = blockstructureinputfunction(Msize)
%
%Specifying the uncertainty block structure
%
disp(' ')
disp(' ')
disp(' ')
disp(' ')
deltastrchkcnt = 0;
deltatoscalar = 0;
scalartodscales = 0;
while deltastrchkcnt < 1
disp(' ')
blkstructure = input('Give the Uncertainty Block Structure in following 2 x 2 matrix form');
[blkrow,blkcol] = size(blkstructure);
        if (((((blkrow == 2) && (blkcol == 2)) && ((blkstructure(1,1) + blkstructure(2,1)) == Msize)) && (blkstructure(1,1) >= 2)) && (blkstructure(2,1) >= 2))
            if ((blkstructure(1,2) == 1) && (blkstructure(2,2) == 0))
                deltatoscalar = blkstructure(1,1);
                scalartodscales = blkstructure(2,1);
                deltastrchkcnt = 1;
            else
                disp(' ')
                disp('Error: The program cannot handle that kind of block structure for now or the block structure is not compatible with size of M(s) ')
                disp(' ')            
            end
        else 
        disp(' ')
        disp('Error: The program cannot handle that kind of block structure for now or the block structure is not compatible with size of M(s) ')
        disp(' ')
        end
end

%return the obtained values

return;