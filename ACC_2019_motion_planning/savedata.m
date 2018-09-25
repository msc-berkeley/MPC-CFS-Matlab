n = size(x);
noo = size(x_c);
x_c(noo+1:noo+n,:,:)=x;
y_c = [y_c ; y];

%%
save([ 'data_2_.mat'],'x_c','y_c')