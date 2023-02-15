
undef("dzdt_lev_merge")
function dzdt_lev_merge(dzdt_infile,ntim,nlat,nlon)
local wwnd_1000, wwnd, dzdt_infile,ntim,nlat,nlon
begin
wwnd_1000=dzdt_infile->DZDT_1000mb
time=wwnd_1000&time		; 40
lat=wwnd_1000&latitude		; 1501
lon=wwnd_1000&longitude		; 2000
lev=(/10,20,30,50,70,100,150,200,250,300,400,500,600,700,850,925,950,1000/)
plev=new(18,double,-999.99)
plev=(/1000.0,2000.0,3000.0,5000.0,7000.0,10000.0,15000.0,20000.0,25000.0,30000.0,40000.0,50000.0,60000.0,70000.0,85000.0,92500.0,95000.0,100000.0/)
wwnd=new((/ntim,18,nlat,nlon/),typeof(wwnd_1000),wwnd_1000@_FillValue)
wwnd!0="time"
wwnd&time=time
wwnd!1="plev"
wwnd&plev=plev
wwnd&plev@units="Pa"
wwnd!2="latitude"
wwnd&latitude=lat
wwnd!3="longitude"
wwnd&longitude=lon

wwnd(:,{1000.0},:,:)=dzdt_infile->DZDT_10mb
wwnd(:,{2000.0},:,:)=dzdt_infile->DZDT_20mb
wwnd(:,{3000.0},:,:)=dzdt_infile->DZDT_30mb
wwnd(:,{5000.0},:,:)=dzdt_infile->DZDT_50mb
wwnd(:,{7000.0},:,:)=dzdt_infile->DZDT_70mb
wwnd(:,{10000.0},:,:)=dzdt_infile->DZDT_100mb
wwnd(:,{15000.0},:,:)=dzdt_infile->DZDT_150mb
wwnd(:,{20000.0},:,:)=dzdt_infile->DZDT_200mb
wwnd(:,{25000.0},:,:)=dzdt_infile->DZDT_250mb
wwnd(:,{30000.0},:,:)=dzdt_infile->DZDT_300mb
wwnd(:,{40000.0},:,:)=dzdt_infile->DZDT_400mb
wwnd(:,{50000.0},:,:)=dzdt_infile->DZDT_500mb
wwnd(:,{60000.0},:,:)=dzdt_infile->DZDT_600mb
wwnd(:,{70000.0},:,:)=dzdt_infile->DZDT_700mb
wwnd(:,{85000.0},:,:)=dzdt_infile->DZDT_850mb
wwnd(:,{92500.0},:,:)=dzdt_infile->DZDT_925mb
wwnd(:,{95000.0},:,:)=dzdt_infile->DZDT_950mb
wwnd(:,{100000.0},:,:)=dzdt_infile->DZDT_1000mb
return(wwnd)
end

undef("locate_min")
function locate_min(data,time)
local track, indices_ctrl, data
begin
indices_ctrl = ind_resolve(minind(ndtooned(data({time},:,:))),dimsizes(data({time},:,:)))
track=new ( (/1,4/), double )
track(0,0)=time
track(0,1)=data&lat(indices_ctrl(0,0))
track(0,2)=data&lon(indices_ctrl(0,1))
track(0,3)=data({track(0,0)},{track(0,1)},{track(0,2)})
return(track)
end

undef("track")
function track(data,time_size)
local track, indices_ctrl, data
begin
track=new((/time_size,4/),double)
do indx = 0, time_size-1,1
time = data&time(indx)
indices_ctrl = ind_resolve(minind(ndtooned(data({time},:,:))),dimsizes(data({time},:,:)))
track(indx,0)=time
track(indx,1)=data&lat(indices_ctrl(0,0))
track(indx,2)=data&lon(indices_ctrl(0,1))
track(indx,3)=data({track(indx,0)},{track(indx,1)},{track(indx,2)})
end do
track!0="time"
track!1="pos"
track&pos=(/0,1,2,3/)
track&time=data&time
return(track)
end

