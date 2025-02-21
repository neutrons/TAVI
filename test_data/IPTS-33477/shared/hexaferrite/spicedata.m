function [data,headertext,headers,defxname,defyname,defxvalue,defyvalue]=spicedata(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function will parse a SPICE data file and return the data in a 2d
% array along with the full header text, the headers in a 1d string array,
% and the default x and y names and there column values
%
%  ML 5-April-2003
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(filename, 'r');
y = 0;
icount=0;
datatext='';
headertext='';
colcounter=0;
nval=1;
ncols=-1;
while feof(fid) == 0
   tline = fgets(fid);
   commenttest=strread(tline,'%s');
   match=strmatch('#',commenttest);
   if match == 1
       lookdefx=strmatch('def_x',commenttest);
       if lookdefx == 2
           defxname=commenttest(4);
       end
       lookdefy=strmatch('def_y',commenttest);
       if lookdefy == 2
           defyname=commenttest(4);
       end
       
       if colcounter == 1
           headers=commenttest(2:size(commenttest));
           colcounter = 0;
       end
       lookcols=strmatch('col_headers',commenttest);
       if lookcols == 2
           colcounter = 1;
       end
       headertext=[headertext,tline];
   else
       tstring=strread(tline,'%f');
              
       % This next set of case structures deals with the possibility that 
       % the number of columns may not always be the same in the data
       % section - this is possible within the scanon structure of SPICE
       
       if ncols == -1
           %This is the firt data line - set ncols to the size of tstring
           ncols=size(tstring,1);
           data(:,nval)=tstring;
       else
           if size(tstring,1) == ncols
               data(:,nval)=tstring;
           else
               if size(tstring,1) > ncols
                   % Coerce the string to be the same size as previous
                   % strings
                   tstring=tstring(1:ncols);
                   data(:,nval)=tstring;
               else
                   ncols=size(tstring,1);
                   % need to now set the number of columns in the data 
                   % array to coincide with the size of tstring
                   data=data(1:ncols,:);
                   data(:,nval)=tstring;
               end
           end
       end
       
       % Uncomment the next line to produce the complete text of the
       % data section
       %
       %datatext=[datatext,tline];
       
       nval=nval+1;
   end
end

%Get indices of defx and defy

defxvalue=strmatch(defxname,headers);
defyvalue=strmatch(defyname,headers);

fclose(fid); 