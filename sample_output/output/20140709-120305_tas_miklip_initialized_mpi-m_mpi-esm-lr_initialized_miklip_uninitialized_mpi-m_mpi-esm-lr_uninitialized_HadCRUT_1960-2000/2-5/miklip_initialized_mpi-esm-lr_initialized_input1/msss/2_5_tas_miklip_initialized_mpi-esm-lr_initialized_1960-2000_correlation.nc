CDF      
      lon    H   lat    $   time       nb2             CDI       GClimate Data Interface version 1.6.1 (http://code.zmaw.de/projects/cdi)    Conventions       CF-1.4                                                                                                                                                                                                                                                             history      LWed Jul 09 12:03:22 2014: cdo timcor /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1960_r1i1p1_196101-196512.ncinitializedZ2LDN2.remapped_ENSMEAN_ANOMALIE_selYear_2-527519BOZIV4Timmean2_5_tas_miklip_initialized_mpi-esm-lr_initialized_1960-2000HX7JDL_merged.nc /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1961-1965Q3A2P4.remapped_ANOMALIE_selYear_2-5277605NLS7LTimmean2_5_tas_miklip_initialized_mpi-esm-lr_initialized_1960-20008RDWS5_merged.nc /tmp/murcss/output/20140709-120305_tas_miklip_initialized_mpi-m_mpi-esm-lr_initialized_miklip_uninitialized_mpi-m_mpi-esm-lr_uninitialized_HadCRUT_1960-2000/2-5/miklip_initialized_mpi-esm-lr_initialized_input1/msss/2_5_tas_miklip_initialized_mpi-esm-lr_initialized_1960-2000_correlation.nc
Wed Jul 09 12:03:22 2014: cdo mergetime /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1961-1965Q3A2P4.remapped_ANOMALIE_selYear_2-5277605NLS7LTimmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1966-1970CH76ZC.remapped_ANOMALIE_selYear_2-527761JUPMVETimmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1971-1975C6WG6O.remapped_ANOMALIE_selYear_2-5277623HN4WDTimmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1976-1980ZGP54Q.remapped_ANOMALIE_selYear_2-5277637Z7LTWTimmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1981-1985T74HT8.remapped_ANOMALIE_selYear_2-527763AIZEZ7Timmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1986-19902WG1Q7.remapped_ANOMALIE_selYear_2-5277601KLFUFTimmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1991-1995IQC7AZ.remapped_ANOMALIE_selYear_2-5277624EBRR2Timmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1996-2000G2YF3T.remapped_ANOMALIE_selYear_2-52776172OAT1Timmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_2001-20053Q7AVD.remapped_ANOMALIE_selYear_2-527763YQWBFXTimmean
Wed Jul 09 12:03:22 2014: cdo timmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_2001-20053Q7AVD.remapped_ANOMALIE_selYear_2-527763YQWBFX /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_2001-20053Q7AVD.remapped_ANOMALIE_selYear_2-527763YQWBFXTimmean
Wed Jul 09 12:03:22 2014: cdo seltimestep,2,3,4,5 /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_2001-20053Q7AVD.remapped_ANOMALIE27763yearmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_2001-20053Q7AVD.remapped_ANOMALIE_selYear_2-527763YQWBFX
Wed Jul 09 12:03:22 2014: cdo yearmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_2001-20053Q7AVD.remapped_ANOMALIE /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_2001-20053Q7AVD.remapped_ANOMALIE27763yearmean
Wed Jul 09 12:03:19 2014: cdo sub /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_2001-20053Q7AVD.remapped /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1986-19902WG1Q7.remapped_2000_obsYT1ZHG /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_2001-20053Q7AVD.remapped_ANOMALIE
Wed Jul 09 12:03:19 2014: cdo ensmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1986-19902WG1Q7.remapped /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1991-1995IQC7AZ.remapped /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1961-1965Q3A2P4.remapped /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1996-2000G2YF3T.remapped /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1966-1970CH76ZC.remapped /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1971-1975C6WG6O.remapped /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1976-1980ZGP54Q.remapped /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1981-1985T74HT8.remapped /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1986-19902WG1Q7.remapped_2000_obsYT1ZHG
Wed Jul 09 12:03:12 2014: cdo remapcon,r72x36 /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1981-1985 /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1981-1985T74HT8.remapped
Wed Jul 09 12:03:06 2014: cdo selyear,1981,1982,1983,1984,1985 /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1981-1985
Wed Jul 09 12:03:06 2014: cdo chvar,temperature_anomaly,tas /home/illing/workspace/murcss/sample_data/HadCRUT.4.2.0.0.median.nc /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc
Updated at 26/06/2014 15:02:28   institution       $Max Planck Institute for Meteorology   institute_id      MPI-M                                                                                                                                                                                                                                                              experiment_id         decs4e1960                                                                                                                                                                                                                                                         model_id      MPI-ESM-LR                                                                                                                                                                                                                                                         forcing       GHG Oz SD Sl Vl LU                                                                                                                                                                                                                                                 parent_experiment_id      N/A                                                                                                                                                                                                                                                                parent_experiment_rip         N/A                                                                                                                                                                                                                                                                branch_time                  contact       miklip-mpi-esm@mpimet.mpg.de                                                                                                                                                                                                                                       comment       Production runs for MiKlip system Baseline1                                                                                                                                                                                                                        
references       �ECHAM6: n/a; JSBACH: Raddatz et al., 2007. Will the tropical land biosphere dominate the climate-carbon cycle feedback during the twenty first century? Climate Dynamics, 29, 565-574, doi 10.1007/s00382-007-0247-8;  MPIOM: Marsland et al., 2003. The Max-Planck-Institute global ocean/sea ice model with orthogonal curvilinear coordinates. Ocean Modelling, 5, 91-127;  HAMOCC: Technical Documentation, http://www.mpimet.mpg.de/fileadmin/models/MPIOM/HAMOCC5.1_TECHNICAL_REPORT.pdf;    initialization_method               physics_version             tracking_id       %4593bd71-f037-4042-93e4-c94e0188a8c8       product       output                                                                                                                                                                                                                                                             
experiment        decs4e1960                                                                                                                                                                                                                                                         	frequency         mon                                                                                                                                                                                                                                                                creation_date         2013-03-17T20:24:04Z                                                                                                                                                                                                                                               
project_id        MiKlip                                                                                                                                                                                                                                                             table_id      Table Amon (12 January 2012) c4004f50cc5c8071abaded9d591b2407                                                                                                                                                                                                      title         MPI-ESM-LR model output prepared for MiKlip decs4e1960                                                                                                                                                                                                             parent_experiment         N/A                                                                                                                                                                                                                                                                modeling_realm        atmos                                                                                                                                                                                                                                                              realization             cmor_version      2.8.3      CDO       HClimate Data Operators version 1.6.1 (http://code.zmaw.de/projects/cdo)          lon                 standard_name         	longitude      	long_name         	longitude      units         degrees_east   axis      X        @  +�   lat                standard_name         latitude   	long_name         latitude   units         degrees_north      axis      Y           .4   time               standard_name         time   bounds        	time_bnds      units         days since 1961-01-01 00:00:00     calendar      proleptic_gregorian         /T   	time_bnds                     units         days since 1961-01-01 00:00:00     calendar      proleptic_gregorian         /\   tas                       standard_name         air_temperature    	long_name         Near-Surface Air Temperature   units         K      
_FillValue        `�x�   missing_value         `�x�   cell_methods      time: mean     history       J2013-03-17T20:24:04Z altered by CMOR: Treated scalar dimension: 'height'.      associated_files      �baseURL: http://cmip-pcmdi.llnl.gov/CMIP5/dataLocation gridspecFile: gridspec_atmos_fx_MPI-ESM-LR_decs4e1960_r0i0p0.nc areacella: areacella_fx_MPI-ESM-LR_decs4e1960_r0i0p0.nc       (�  /l        @      @$      @.      @4      @9      @>      @A�     @D      @F�     @I      @K�     @N      @P@     @Q�     @R�     @T      @U@     @V�     @W�     @Y      @Z@     @[�     @\�     @^      @_@     @`@     @`�     @a�     @b      @b�     @c`     @d      @d�     @e@     @e�     @f�     @g      @g�     @h`     @i      @i�     @j@     @j�     @k�     @l      @l�     @m`     @n      @n�     @o@     @o�     @p@     @p�     @p�     @q0     @q�     @q�     @r      @rp     @r�     @s     @s`     @s�     @t      @tP     @t�     @t�     @u@     @u�     @u�     @v0     �U�     �T�     �S`     �R      �P�     �O@     �L�     �J@     �G�     �E@     �B�     �@@     �;�     �6�     �1�     �)      �      �      @      @      @)      @1�     @6�     @;�     @@@     @B�     @E@     @G�     @J@     @L�     @O@     @P�     @R      @S`     @T�     @U�     @�	     @v�     @��     >��>�A.`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�?��?��`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�l�=f�=�U�>Tȯ�7Q�W5�>���?7<?ʦ>�D=�|�?)��`�x�`�x�`�x�`�x�>��+>�\�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�>���>�%�>�P���H=�����U�>h~�>$��?!�X`�x�`�x�`�x��־`�x�?� ?>`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�?]bc=�>�^���I���b=�G�>k+�?-P?3Wt?��>���>���}Ms�t5K�ĲE���?"^r>�-�>���?/�ҿ
8��cA�bgF:�f+=��>S�`�x�`�x�`�x�?�  `�x쾱���UL+=_��>-<>�Ь>�Nx>�:>����������>���>��>�,?۫?Cb�?B:�>�%�>�0�����gI��%?"f?%������?;�
?��<��?��?�|�L�?�?��?"��>�:_?+�>B�揞?L��?>Q��뱶�����r���m_>�~�=O��>=)�>���?�;?CY�?fф?ef�?O�?7�?GP�?]I>�u�>��G>Mi��p�8;�zs�x�;��m��9S>B��>*�>9�>�e�>��t>�G�>��U=���?Rxa?=?Q�>�g>�l�>"s%>	B��=*G>H�;�} �ثǾ�s�
� ��/� <�L.^���$->�N�>�f>��>���=�*0=�����F�=��p>�0�>�O?�>��=�w�>��~>�V>"Y�<Y�?��?#�u?e^>�K�>�K>n,�>��ʾ�Q���/��n=>HX>oS>�܆>���>�yk>��J�'�>_ *>�:�?|�>�2j>͛>������/�������=D}�>��`<ǩ:?���`��=���=�i
?+�??��?>�?g  ?\1? �>��4�ںS�<�)�?�Lt����l�V$�4�w�6���� ��ཨ�~��ë����=��>�غ+:o����˚
�������,��>N-�?��>��(?	�? �>[rJ>�S�>��L>�)0>�$>��>�>�rr>I�9>���??�?
k�>��#C�>f�����ƿ�̿$��� ���j���-�>K�X=���>�UO��羑��݉n=�JO>�\9<���<Q�>��|>��/>��B? !?A
�?7�
?��>��Ľy%>2�U>��[?	��>��[>|��OH�`ZP��3پ��1��b�K��ph��b������S�<�0�>��F=�S=���=�~==�I�Yq=��>��>�W�>P6:���=�х>���>h�P��5���4=�>���?w�?mR?1}>@�=�Y�>�0>�=>�-�>���>�/]=�A���> �@=0�>ΫB����=^��������Ce>7zK���~�>
?���>$��U��ʑ=z>ii> Y{>�C>ɋ�>�&%>�R�? ��?�-?��?=ɧ?��>�P>��C>�x~=ʓ�>j�J=:w=���	I����"�Պ��m�i���<��E>�ˣ=c��>���?B�?;]<?;��?\&�?2R?�J>�j�>r�>�g�>_">��O?Hm�?	�*��侀2��n�!��>>"�>_�u;Yc:g�����Ǿ.l�<-�>�x%>�i->��f>]3%= ~W?5ݼ>�b�=���/S6=��v>Û�>��>�Gm>ce�>��>믽@,w=�i>� ���о�!�<ގ�>�"�>ŕ�>��>���>���?��>��?��>�
�>]o�>�2�>�P>H�>�:�>�h�>4'y<ꑪ���پ־����߾Wm�=~��Nj����?D^>�t>��?!?uu?$��?\�?NQ?`��?g�D?.�L>�.?' ??L�6?!K�>�P�>��>�N?�<?��?5B�?6�>�{�>��}>�n�?��>��>�U>���=�Z<����h�=��{>k}Q>�.4? �?A�>���>�D==�B�� �o=�t����������ڿ%���/��" >���?#_b?�O>�;>�1�?���w��K>���>��? �>�Q>D�4>�*>m*v>>�>M�=��U�
i?#�? �v?Z�?E�$?&r>�	+>��+?�x?v%? l$?\~/?\�P?@4??��?Ep]?:��?.�?.�?<�?Y�?A�(?-��>���>y>��B?!�>�Z>�eS?6�?5�B?N  ?�'>�V�>��>�e�>���>|��>W8>Z��>�m�>��?�?/>�ǰ>枅>�@s=�h��	i-���3��T4��E���-=y��>J�C��f\>�u�~��>���?!�n?��>�f�>���?F�>��>/����k=f��>Y];���>,�7>��
?	|�?DW�?d1�?c ?>�f?A�P?J^�?[��?iI?N�9?QF?L�??�? F=>�T>�7>x]?
�/?1{�?X�?Dn?4��?6tp>�A?	?AX?*Ӝ>��>�Q�?L�K?Y :?da>?1\�?�)>�U�>��k>}�=�]>��?-?&!X?��?
�0>�ʞ>������L>�:|<}#�
�;��p��=`���b��f`>���?$:>�.>�{=�s�?))>R��>���?)��? �>B��=���>�Y>x�>MG">��>�&?�*?-��?,6?WO�?m��?P�j?+B�?b}?D�)?g)?q}?\�?J��?�>ظL>z>>lI~?�K?`�w?r��?l��?B�>�S�>�t>��*??	�c>Ʌ�>���?�`?7��?8��?P��?&>�^?&1�?m�>���?��?6��?R�?��?��?YE>��>���>��x>�o>ђX>H�>>���?'$? ��>��2?�>��2=��i��J=�@'B���#>�����1k:=�@�?
�?W��>�cA?@U?W{?o��?t/D?ng�?a-�?;?�?C�?@V?FsA>�#?>��?	��>���?=�?Us?R��?6��?&0r?-�???3�?B�?Q�y?b�i?JŞ?.��?8�?	�I>ѾX>Ղ.?>��?(�>��?�?�>���>���?
W?��?>?%�i?$v�?�?)��?��?�'>���>ގ�?�f>�G?�?%/>�J]?�U?N��?D�N?�?�g?�>���>�S��wLp���.�1��?`>�c��-���>��j?B�W?ItZ?BSK?T��?:k?��?5�>4��>g��>��q?;��?'*>87-�c��?^7Y?N+�?=�?�?:?�?;\�?MI�?Ak�?Q'
?Y��?V��?,(?�i>�4<?�w?^>>��2? �"?(ɸ?2j�?=��?X�>�χ>�+>�c>�܉?��?�?/K]?C��>���>� ?(r�?(�3?��>��?	��?I:?/Yh?3'?�_>�j�>��E>�_�>ߓ�>�8>���>'\�=�i�����/���F�'�=L)�>c�������U�?�?c��?Xl�?@��>ԯ�<Vc>1�> �׽%	g>��?&��?�j?!1�>��?�$?{��?]\.?Q��?!.�?<�t?cd�?J��?$#?(�+?0<?0�~?#7�?%�?)�R?*u2?1k�?=��?]nj?g�?df?Q�?3�?+�?W�?C>�0R>��?@Z�?D�)?,�r?$%�?1q�>魶>��>�XA>�>�T�>�(�>��? T�?��>��B>Ϩ >��@>i��>X��>$3{=��l�}TE�f��G}�&ؾI�E>9"E��f�O˫>��V?%�0?2/?;9�?N3?7)>~��<���,Q�>��? \�?�?_�?��>�,�?@�q<.|_�d�?U�?3��?N�@?P�B?W�h?U�?IR?G�s?B>�?A
?;��?9y�?/�]?OU�?h�l?[�?hvr?v�}?m��?f3?bw�?dyv?K��?!��?F�?/�?D��?2�%?�>��>�9�>�G�>A��=��>h�>���>k�>4[>B�>���>��>�U�>�=.Q��*ˋ��--���#o�<l�G��=籠>
��m�Ҿ�Ym�̟Q`�x�>���>��>�_�?I4�?0�n>N8�>OIJ?}?.>���?@+	?;��?/&�?[p�?I��?N�?J(�?J
�?a��?Z6?W�$?OY�?D��?Mu-?Xe�?]S9?V��?Q"?F$�?L�?*�Z?]mR?w�?eS�?I��?w�@?g�!?_IJ?��?%�?b]>�S�>:�>s��>�o$>~�c><6a>%�A=oL����x�	P��i���i�7a�{?��V�1B�j;�&��Öd��������?���ܴ�ƌH�pi�>�#���؟�Dz>���>��d?{�?1ֈ?*�<?%�1?Ew�?R�K?@:�?7ǹ?=�?89f?0��?x�?oj?m�?v'?LH�?8:�?d��?V��>Gu�>���?X��?�?"��?Q��?W�?_�>?`a?c#G?_��?_H�?Up-?iB�?[{�?V�Y?,I#?g��?Rڧ?'`�?.
�??�?��>��l>�%=>�=�Fb>/��>B��=وJ<,)D����A�}��;�AJ��$����|�o>l=��f>��=@�N=�a^�&t�=.E�>;k6>5b�>Q͵>D�X=�I�>�?55;?	f�?1�?�?0Z�?�:?*�?N�d?gq�?h[R?Ws�?M�N?<�6?2�i?@*? �(=?hi>�$>��K?Ev#?W��?D�?R	K?a~�?5+?P�M?>8�?"�k?Z5r?Y��?V��?Tu�?B�|?M�?I�h?c�?a��?d��?_��?h �?r|�?n��?r�$?k��?eU�?!ɾ=���>\�>�� >Pp7>��i>�� >A��=�L�i�<���~�x=��=�_N=���>c��>��`>��V>��d>�ˮ>�}->�z>�,x>�F�>��5?�?*�U>�Dj?>7?F`? �p?�O?4W�?T��?V�?\�.?e�.?o0�?ZZT?\ �?J_�?jO7?Q��?J�i?=�d?b=?]��?L�{?O�?_u�?M�?LF,?:(C?7B%?�%?�b?F�O?Ew?M�T?)�?T�2?_
s?`��?K,m?RK�?Vf?`�?jH/?%��?Bb_?d��?Q8�?��>��1=�=�=y$%=���=$q����v�R�V�4^I������8]5��3�<���;��L��=�=�6�>U�d>�y~>��
>�bH>Lv">��?�*?�F?J*)?K�??�b?L��?DFG?.f?+B?7o�?DH?WN+?e	�?a�?tӔ?\�I?k��?l�?]��?Z��?V� >�ۿ���j`�x�?e��?q�g?7�=>���?#�?3%�?:�,?(�n?#L0?9c?3nP?Y�?>>?j�?�?e��?`^P?3�0?�J?8��?X?Xq&?Q�?k�?AW>��M>C�E>KV�=� �>e>���>r�t>��<>3><=-�]=�[m>+0G>�%�����'V������~��`-�};�>�P�?+3?0��?o@?�?�d?3�v?4%?)��?;j?:�7?!�?%��?8�?N�@?U?N�B?>?B4�?q��?\�?b1�?kHt?Y<�?,�o`�x�?'io?5j?L�?[�?({h?D5�?IP�?3Mf?1L�?<�z?A�?S �?SA�?0��?ǒ?Jh?7��?@��?FB�?a�@?G�T?o�?V??OJ!?f+@?js!?j�?@��?��?#��?.T?�P>��>���>,r>,Z&=�`'�vd�bh�=·-=��콁��w w��;U� �,�Os>��?P�?cJ%?K��?@� ?Dp\?2�N?��?�x?8�?G�I?1Za?��?<�?O�??�N?CN?L�A?=�9?t��?f?`�x�?^��?D�n?V>6>�&�?e\>�8G?#ƥ?+D?X�?BT~?�)?C��?Pf?KM�?`��?<?�$?&-�?W��?_M?q�u?Qŗ?8��?U¿?_ʊ?Ks$?:��?QE�?V��?j�&?h�x?S� ?S(�?P޸?>(h?*o,?��?&��?��>�t>��=��\���L�4�(:ž\�S>4q��&ʸ����=U͂>ە�>�XZ?W�J?q�@?e��?A�?#M�?-?6�?F#?B�g?Wqh?O�#??hd?>�U?K��?F��?W��?g�?b��?v�?_6T?G�?TG3?X�?O�?0��?�]>鬸?��?74A??�V?0C�?Z�^?S�B?I�b?H�>p+V>2�@>�!�??�`?"��?a�?Y�f?Z�I?d��?sA�?f�?Z�Q?g�?^�(?K4g?Z3?LO??�O?OEY?NE�?I�<?Q&??�^?*�f?#�-?
�>��p>�g�>b�n�q�Ⱦ�J<��<�����?.7�?Q"E?F�?+#�?M��?c�?E��?!>�Ce?��?e�?��?R:�?[ޠ?Qݡ?L��?>�^?$��??��?@�?[V?\�d?VO?e5�?YE�?\��?X�i??s>6
F>u�>���>�P?O�?/��?W�f?k��?^�?F�?9O�?J��?T�"?q��?v��?_Y*?DX�?`J�?c)�?]�a?h�?B�?H�m?$�#?f,>��>�ch>[�;�m���z�>D��>�đ>���?�?�3?ȡ>�>���>�>�>��>�[�>�,�>�5?<��?P$K?G��?4I�?"�?)��?=tW?<�O?#Ç>��=?m?C�w?Ar�?Uw�>��J?��?R�&?r ?8�?%u�?L��?q�<?W�/?P�L?l�e?xj�?m�a?Qf?K'>?W���c��>z|�?�?B?^�h?o�:?h$�?^�&?VV|?I&�?G�??1��?_9?^p?dS?]��?N��?<+�?- ?�?5y�?
�>=�:�I������ؾ���� Š��ݾ�>���>��>���=K�s��'��}i>أ>��y?K�?)�?8��?A�?(�?6��?)��?�?��?#�?,,K?1��?9�w?��?Pa?>$N?8�#>Ph,?
�>���>�(�>�g?
�:?Z �?a��?]�?Y�}?>M5?Ro�?V��?.��?
��>�i�>�ֱ>)�>�R?�?'Y�?��?4�w?Dh�?$�?P�?��?��?C��?Z�P?\��?JP�?=,r?#� ? ?q�?^�??��?Br?G�?1�$?#�x>�<x�>�=�6�=}D�=��$��*s�������W	>^3?8�.?]u�?[�%?K�m?T��?;da?0?Y?&J�?��>���?��?#�J?1�~?8�,?(�<?)b??�&?8�k?6�>��2=�yt>�#�?�>��?H�?fɀ?P�3?Sf�?@->?߰?��?�J? �>��(>�0%>���>�K�>�(�?z�?J�?*?6��??/�?0��?(\A>�8g? ar?6!�?K�??J[�?F'?@P?M?�.?��?28�?BJ?Os�?dzO?6��?*��?��>��>J�+=���<��:�4�T��;f1#=��>T'>>�/�?�?+�?$D?;I+?>��?0�M?&�?�?��?��?D�>��z>��d>�E?s
�>�@�>���>��]>���<�芽����G��>Aʡ?{�?3�8?(8�?<�?�?H�>�EF>��3>g�>��c>���>�
>�n>���? h�?�t?�4?�.?��?�t?.�{?)��?!�?��?�?<O?T�?�v>�5?� ?:}?;�?M�:?H3�?3��?&HU?-72?-^J>Ћu>�=o1K��8;��>)&�=�hW>6RT>�R�>�T�=�+
>.5.?�Z?B�q?Mƒ?B<�?8�T?'�/?N�?�?	'�>􆼽M�a>�}�>��>�'�>�b>��/>HR�>�\�N��1G�����]��+lžH��=��>�yt?�>�N<>�5>�b�>�l>��>�(`>���>�q�>c�E>d��>Y��>|m�>�<�>���>�C�?e\?|?9�j?+>?��?�>��>�9>��(?	R�?:�?^-W?U|�?@��?;�??PZ?8��>��m>F�A�c�ᾍmI��.Y�7G�>77$>�V�>��>��1?��?"��?,�M?22�?)�z?9&�?@��?6O�?�K?y�>��t?��??ک?�E?,>>�1�>'�>��>��>�7����<� e>Z
�0�辏k����I��]�>8�R>�KB>�� �-�)=�.�=q.>>�.$>��R>+V�<�ee<�#=ԏ�<:��<b^�>'��>�=v>��b>���>�
/?�.?.{�?<Z?7(�?P? D>�{�>���?x?$�?J~�?*�Y?.D�?8�e?5�?��>�x/>�;4�&<�"}=}M�=]Y�>��1>�j?J??��?4�3?��?/�H?H�?O^�?@�b?*�?*��?$��?"��>6/?>>M�??�S?6CN?k>�a�>� �>�Z>>ku�=H����v`�x�>���>B H�Bx�=���?;V?5�0>�v�=_��>�n>#��>��D>�w�>Ly�<�]��j�r���w��t��� ��yb;	�X>\�>�f�>��`>��S?J�>�>��>���>�e	>��
>td>Vً>�ƻ>�{�>���>���>k5X>/�>X!�>��>�7*>��X>���>Қ�>�)�� 9J>��?;B!?>~�?=�?S�?N��>�of?-{j?1<�?1�`?Ah#?Ҳ>�T�>�s4???7P�?N?
�X?M	�>ϣ>ɘ=�u�=v�>`�i`�x�`�x�`�x�`�x�?T�?Es��4?I�f?PPn?�7>컫>�͞>�2�=�"ܽi�:��<־��f��K�����=�t>�c�>Q��>p�Y>��I>g<>�#>Ǳ�>� �?ђ>�Z>���>���>�c��Dq;�ӏ���>�t�>�'�`�x�=��Ѻ�ݩ`�x쿀  ��  �%㾹�;��  �>.��?:>F���-��.��  ?�  `�x�`�x�`�x�?Sm)?U
��  >���?����]��Ž�]�>g}k>��>��>�|�>d8�>�7m`�x�`�x�`�x�`�x�`�x�`�x�?Z>�?I[˾�2T�"��>���ݟ>V�۽����ގ>�H ?4Ł�E�VCſ+�׾YZ�>� �>���>���?e�? �>���mv�\Z��Z� >��8?�  `�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x쿀  ���Ŀ�  `�x�`�x쿀  `�x�`�x�`�x쿀  ����ٿ�  ?�  `�x�`�x�?�  >���`�x쿀  ��  >��<?V�;?m?
/>
	$`�x�`�x�?	��?X�`�x�`�x�`�x�`�x�`�x�`�x�`�x�?^�S?]��?.g]>�S�`�x�`�x�`�x�`�x�`�x쿀  `�x쿀  �Sȿg���x㴿y/�w�&��  `�x쿀  `�x쿀  `�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x쿀  ?�  ?�  `�x�`�x쿀  ��  ��  `�x�`�x�`�x�`�x쿀  ��  ��  ��  `�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�?�  ?�  ��  `�x�`�x�`�x�`�x�`�x�`�x�`�x�`�x�