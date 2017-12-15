function imagetoprint=printLocal(image,filename)
   
   % if nargin<=2
   %     viewimageFolder=pwd;
   % end
   % currentFolder=pwd;
   % cd(viewimageFolder);
   % image=viewimage(image);
   % cd(currentFolder);
   % close;

   %Print EPS
   figure,imshow(image,'Border','tight','InitialMagnification',100);
   filename_EPS=[filename,'.eps'];
   print(filename_EPS,'-depsc2','-r300');

   %Print TIFF
   filename_TIFF=[filename,'.tif'];
   option_plot=get(gcf,'PaperPositionMode');
   set(gcf,'PaperPositionMode','auto');
   print(filename_TIFF,'-dtiffn','-r0');
   set(gcf,'PaperPositionMode',option_plot);
   
   close;
end
