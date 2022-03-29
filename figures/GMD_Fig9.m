load mass_balance_inlet


for ii=1:length(param1),
    

    for jj=1:length(param2),
        out = output{ii,jj};
        Vbarrier = b_struct.dy.*(sum(out.d_sf.*out.d_sf./out.s_sf_save./2.*(1-(b_struct.s_background./out.s_sf_save)),1)+sum(out.h_b_save.*double(out.x_b_save-out.x_s_save),1));
        
        blub(ii,jj) = Vbarrier(end);
       
        Qinlet(ii,jj) = mean(out.Qinlet);
        Qoverwash(ii,jj) = mean(out.Qoverwash);
 %      plot(out.Qoverwash)
    end
end

Qoverwash(:,1) = Qoverwash(:,1)*3; %correction for doing part of the run
F = Qinlet./(Qinlet+Qoverwash);

%figure
%histogram(F(:),[0:0.02:1])
%
% F as a function of grid/timestep
c = {Qinlet./Qinlet(6,4),Qoverwash./Qoverwash(6,4),F./F(6,4),bsxfun(@rdivide,blub,mean(blub))};
ct = {'Qinlet','Qoverwash','F','Mass'};
for ii=1:4,
figure, hold on
title(ct{ii})
imagesc(c{ii},[0.8 1.2])
colormap(jet)
[a,b] = meshgrid(1:length(param1),1:length(param2));
text(a(:),b(:),cellfun(@(x) num2str(x,3),(num2cell(c{ii}(:))),'UniformOutput',false))
set(gca,'XTick',1:length(param1),'YTick',1:length(param2))
set(gca,'XTickLabel',cellstr(num2str(param1')),'YTicklabel',cellstr(num2str(param2')))
end
