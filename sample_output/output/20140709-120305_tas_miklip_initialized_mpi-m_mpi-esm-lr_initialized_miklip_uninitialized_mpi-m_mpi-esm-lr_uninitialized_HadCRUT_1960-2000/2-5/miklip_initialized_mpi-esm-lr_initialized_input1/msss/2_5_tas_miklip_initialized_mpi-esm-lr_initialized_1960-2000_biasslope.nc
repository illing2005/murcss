CDF      
      lon    H   lat    $   time             CDI       GClimate Data Interface version 1.6.1 (http://code.zmaw.de/projects/cdi)    Conventions       CF-1.0     history      �Wed Jul 09 12:03:23 2014: cdo mul /tmp/murcss/output/20140709-120305_tas_miklip_initialized_mpi-m_mpi-esm-lr_initialized_miklip_uninitialized_mpi-m_mpi-esm-lr_uninitialized_HadCRUT_1960-2000/2-5/miklip_initialized_mpi-esm-lr_initialized_input1/msss/2_5_tas_miklip_initialized_mpi-esm-lr_initialized_1960-2000_std_ratio.nc /tmp/murcss/output/20140709-120305_tas_miklip_initialized_mpi-m_mpi-esm-lr_initialized_miklip_uninitialized_mpi-m_mpi-esm-lr_uninitialized_HadCRUT_1960-2000/2-5/miklip_initialized_mpi-esm-lr_initialized_input1/msss/2_5_tas_miklip_initialized_mpi-esm-lr_initialized_1960-2000_correlation.nc /tmp/murcss/output/20140709-120305_tas_miklip_initialized_mpi-m_mpi-esm-lr_initialized_miklip_uninitialized_mpi-m_mpi-esm-lr_uninitialized_HadCRUT_1960-2000/2-5/miklip_initialized_mpi-esm-lr_initialized_input1/msss/2_5_tas_miklip_initialized_mpi-esm-lr_initialized_1960-2000_biasslope.nc
Wed Jul 09 12:03:22 2014: cdo timcor /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1960_r1i1p1_196101-196512.ncinitializedZ2LDN2.remapped_ENSMEAN_ANOMALIE_selYear_2-527519BOZIV4Timmean2_5_tas_miklip_initialized_mpi-esm-lr_initialized_1960-2000HX7JDL_merged.nc /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1961-1965Q3A2P4.remapped_ANOMALIE_selYear_2-5277605NLS7LTimmean2_5_tas_miklip_initialized_mpi-esm-lr_initialized_1960-20008RDWS5_merged.nc /tmp/murcss/output/20140709-120305_tas_miklip_initialized_mpi-m_mpi-esm-lr_initialized_miklip_uninitialized_mpi-m_mpi-esm-lr_uninitialized_HadCRUT_1960-2000/2-5/miklip_initialized_mpi-esm-lr_initialized_input1/msss/2_5_tas_miklip_initialized_mpi-esm-lr_initialized_1960-2000_correlation.nc
Wed Jul 09 12:03:22 2014: cdo mergetime /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1961-1965Q3A2P4.remapped_ANOMALIE_selYear_2-5277605NLS7LTimmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1966-1970CH76ZC.remapped_ANOMALIE_selYear_2-527761JUPMVETimmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1971-1975C6WG6O.remapped_ANOMALIE_selYear_2-5277623HN4WDTimmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1976-1980ZGP54Q.remapped_ANOMALIE_selYear_2-5277637Z7LTWTimmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1981-1985T74HT8.remapped_ANOMALIE_selYear_2-527763AIZEZ7Timmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1986-19902WG1Q7.remapped_ANOMALIE_selYear_2-5277601KLFUFTimmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1991-1995IQC7AZ.remapped_ANOMALIE_selYear_2-5277624EBRR2Timmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1996-2000G2YF3T.remapped_ANOMALIE_selYear_2-52776172OAT1Timmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_2001-20053Q7AVD.remapped_ANOMALIE_selYear_2-527763YQWBFXTimmean
Wed Jul 09 12:03:22 2014: cdo timmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_2001-20053Q7AVD.remapped_ANOMALIE_selYear_2-527763YQWBFX /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_2001-20053Q7AVD.remapped_ANOMALIE_selYear_2-527763YQWBFXTimmean
Wed Jul 09 12:03:22 2014: cdo seltimestep,2,3,4,5 /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_2001-20053Q7AVD.remapped_ANOMALIE27763yearmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_2001-20053Q7AVD.remapped_ANOMALIE_selYear_2-527763YQWBFX
Wed Jul 09 12:03:22 2014: cdo yearmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_2001-20053Q7AVD.remapped_ANOMALIE /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_2001-20053Q7AVD.remapped_ANOMALIE27763yearmean
Wed Jul 09 12:03:19 2014: cdo sub /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_2001-20053Q7AVD.remapped /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1986-19902WG1Q7.remapped_2000_obsYT1ZHG /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_2001-20053Q7AVD.remapped_ANOMALIE
Wed Jul 09 12:03:19 2014: cdo ensmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1986-19902WG1Q7.remapped /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1991-1995IQC7AZ.remapped /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1961-1965Q3A2P4.remapped /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1996-2000G2YF3T.remapped /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1966-1970CH76ZC.remapped /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1971-1975C6WG6O.remapped /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1976-1980ZGP54Q.remapped /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1981-1985T74HT8.remapped /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1986-19902WG1Q7.remapped_2000_obsYT1ZHG
Wed Jul 09 12:03:12 2014: cdo remapcon,r72x36 /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1981-1985 /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1981-1985T74HT8.remapped
Wed Jul 09 12:03:06 2014: cdo selyear,1981,1982,1983,1984,1985 /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1981-1985
Wed Jul 09 12:03:06 2014: cdo chvar,temperature_anomaly,tas /home/illing/workspace/murcss/sample_data/HadCRUT.4.2.0.0.median.nc /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc
Updated at 26/06/2014 15:02:28   source        CRUTEM.4.2.0.0, HadSST.3.1.0.0     institution       LMet Office Hadley Centre / Climatic Research Unit, University of East Anglia   title         AHadCRUT4 near-surface temperature ensemble data - ensemble median      	reference         �Morice, C. P., J. J. Kennedy, N. A. Rayner, and P. D. Jones (2012), Quantifying uncertainties in global and regional temperature change using an ensemble of observational estimates: The HadCRUT4 dataset, J. Geophys. Res., doi:10.1029/2011JD017187     version       HadCRUT.4.2.0.0    ensemble_members       d     ensemble_member_index                CDO       HClimate Data Operators version 1.6.1 (http://code.zmaw.de/projects/cdo)          lon                 standard_name         	longitude      	long_name         	longitude      units         degrees_east   axis      X        @  0   lat                standard_name         latitude   	long_name         latitude   units         degrees_north      axis      Y           p   time               standard_name         time   	long_name         time   units         days since 1850-01-01 00:00:00     calendar      standard        �   tas                       	long_name          near_surface_temperature_anomaly   units         K      
_FillValue        �I��   missing_value         �I��   reference_period      ��     (�  �        @      @$      @.      @4      @9      @>      @A�     @D      @F�     @I      @K�     @N      @P@     @Q�     @R�     @T      @U@     @V�     @W�     @Y      @Z@     @[�     @\�     @^      @_@     @`@     @`�     @a�     @b      @b�     @c`     @d      @d�     @e@     @e�     @f�     @g      @g�     @h`     @i      @i�     @j@     @j�     @k�     @l      @l�     @m`     @n      @n�     @o@     @o�     @p@     @p�     @p�     @q0     @q�     @q�     @r      @rp     @r�     @s     @s`     @s�     @t      @tP     @t�     @t�     @u@     @u�     @u�     @v0     �U�     �T�     �S`     �R      �P�     �O@     �L�     �J@     �G�     �E@     �B�     �@@     �;�     �6�     �1�     �)      �      �      @      @      @)      @1�     @6�     @;�     @@@     @B�     @E@     @G�     @J@     @L�     @O@     @P�     @R      @S`     @T�     @U�     @��P    >��->ð�I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I��?�{F?�G��I���I���I���I���I���I���I���I���I�ʿ���=���=��>QK#�_�4#>��g?E�? �g>H��=�[>���I���I���I���I��?4K>��*�I���I���I���I���I���I���I���I���I���I��?'��?��>�f��C�=�/�͔>�j>���?��9�I���I���I�ʿxb��I��?6b�?<�@�I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I��@X�`>'�Y>�ֽ�����?<��o=�G�>�*?B�?;G�?��>�~r�����؎�7^�`�E>�"�=��>;*7>�R�x,�X{��(�;d�=~��>���I���I���I��@^3[�I�ʾ���/O%=

�=�I
>94W>�
�>�(>�U;�3�\�A_?�>x>�6a>!j?Fy�?~R?qq1?'g?�^��������"?�m�?n"�u??�]_?�yR=BS�?�?�?��}��Zk?*�6?�@@?��4?���@$?K������?�$>��`���ɿ�T��㾄�w>aN)<ڄ=��>O��>���?b�?d�?=]7?��>��L?Y1�?W��?�7>�V9>������;'Ñ�
ꧾ�C�S��=�r�=��>"nz>�)?�u?4�?cb=�E^?��?�3?Vp>�݆>���>PN>4z���<=��X>�j<,u��4���6¿WH�����ے��i��@fݾ9�N>���>��>��*>��=��=��c��u�=�|>�W�?DK�?*�>�h1>�?9�?tW>7P4<��t?A�->��?:Ȗ>�3s>�ڹ>�&,>��v	�x�Ǿ�����#�>��>��?jy?�?6P�?}��j>=��>��?l�?a�>�\�>�-�� �����u����=1;�>h_-=Ǟ?R�0��[Z=��->�C?E��?Hʝ?f�3?��?|T�>��>�eؿѿp�ѿ�̽����T��������!Pȿ�hq��д�?	��s�&��>mk�>��ϻ�i���E�|ö�{T��j����>��?��?D��?���?��q>jY�?��>�>��?O�'?���?Z�p?��>͆ ?�q?�<?���?�?�T��j>�r3��������\�^%��)��h3/>���>$�>�3N���C���v��>6�~?��=d�<��4>���??�?N�?Ue�?���?�C~?���?�;��^>�՘>��>�T�=��>
��Gk˽�����¾��ܿ����Wh��ep�mꏿ�M�֬�=�?5W�>Xf7=�c>8�>	&+���8>>�f>��4>��X>sIZ;>hA>���?�N�?.����%Qp=�:	??9�?�/.?���?���?G>�_s?�@�?���?P�?m1�>��>Ӡ�5bx>��G=���?O�S�.�=ʑ>�6%<�QMf>�Y뾂v�����^J�>Z6�}�C�:/o>N}>��[>�h�?~;1?���?\)�?x(?M�K?>)?3��?�)Z?���?T0?�x>���=�t>�2>�2>,��������C�وQ���/���eJ�=`\�?=\2>��?eo*?�l?�mA?��?���?�1�?���?+��>�&>��i��$?�o�@o�'?:�e�N�־�������AA�>���?	z�;���:�J+�,���t�<B_=>�K�?>1�?I�l>��<�$�>�:�>���=���$�>�)?��?<5g?�@>��*>�-�>0$F�e�?=��>��V���!��S�>��?�@
!?�p.?�a�?#'�?+2�>Q��?F{�?�L ?��@�^?��h?& M?T�?%�>�?=����:����dc���>E?����*��?�Q�?_��?`�?�:&?��|?���?��j?�.�?�Q�?�Q?�"�?q��?i?yש?a��?�Q?ɨ?Rڏ?2q�>�1L?[}?B*~?�B>z��>� �?#U�?N�t?1��>��>
-�<<޾׾=�k�>۳?g�?X�?�>Ӎ�>�=�U��1��><�9��l���v������-���_?7�2?Þ,?�j�?4R�>��j?y���Ͽ�c?��?��M?��=?7h>���??��>�zE>��>��
>@վ%�?�4?�z?x�?�wk?.,�>�&'>�Jv?+b$?A�?Yf�?�c�?�Z�?��V?��%?�K?�R�?�y?�Q6?͛?��R?Œ�?�-�>�e�>%��>��>�
�>���>{F>��>�7r?\�3?Nc�?*�a?V? �>�+>�"�>�]^>��?)D�?<�3?S)�?=(�?/�F?P{?(^�>  ��4��|9�n�ݿ�ye�<Ν=�LQ>l�u��\�>%�Ⱦb��>� �>���?�->y��?C�?��@>��>t��*��=專>[�L;>S��>�.�>�BT?�&�?���?�p�?�$k?9B�?�2?��?��^?�/�?}�M?x��?�Z?e�?�>��>©b?[\/?�E�?�|{?�B"?Ѳ�?\�>��>�?�>�|U>��>L�>EmQ>ӈ/?�o>�	9?Y�?X9f?�4?@Tm>��R>:5z?y�?��?��$?|�6?qo�?ff?#5�� ����?�C<�1��\��3�˾:&���*r��Ѣ>g��>�}!>�`�?@X>��>�jj>,�>�v�>�m�??�">�k�>?�>ݔg?9�?�?{��?�k�?�
�?d�??�?O*�?R�G?I >��>�c?e?���?��_?��x?v��?q*,?1��>^_�>�j�?lD?�;�?�U�?��]?��>�u�>i*�>�w>�2�>�,T>d�<>�T�>��?w�>�<e>�M�>�fS>�i�?6��?G��?�R?L`|?�3?e`6?Xf?>��?Id?Ipc?4�j??p%?PkQ?�>�J+>�?K<�>�	>�ד>��m>3kt=8^���~S�1E���{�>`k�����1=��l?3 �?�O�?�I?r#�?� ;?��0?���?��~?���?d��? yK?�#>�#>���?�#>Ρ>�f�>�&�?7�?9Ҁ?e�?bC?�0�?���?E�Q?-@�?*[ >�,�>�^>�F�>�p?>��r>�&�>�t">�٘>�-�>���>���>I|�>I<�>��>݂�>�s�?:�?D%�?H�E?�;�?K�!?+�v?�s?p�?T�?>٦?#�~?Bq�>���?D��?���?[��>�BB>߉4>���>�+?>*Ľ'��J;-�&�����>��"��9����>��?\2�?5�`?>*P?p@�?Q\d?*��?� >-�$>�y@>�X�?
W�>ւ=������?X�
?k7?��>�Nq?�>��L?�?%g�?H�1?;��?��>�\�>���>rJ>�[�>�E�>��z>�?W>�Y�>�P�>d�>,�R>TC�>�ʥ>�	>��&>��\?C��?n��>䑭>�j�?#z@?C��?)d�>��? ?��?-a;?L��?$�Q>��>ʁ�>���>˿>�ٔ>���=�c=�t��y}U��[�����z=��=�i6�t�D����>�=�?h<?g^�?2X>���<1�5> -`>$��.�>u�??p�>��>�GB>�Jx>�+�?R�Q?��?V1�?I6�?]h�>Ѥ�>ߵ�>��>���?S?��?�?|?5N`?.Ju?��?!�?'z�?0ט?��>�;�>��>��?��?�p>�bA?��?Iқ?3�n>� r>ì<>�v�>�>�RV>�ee>��>��>���>��
>�۱? T?>�[Z>�P@>Z��>@��>$ �=��?������K��4�U�@������>0ƻ�yQ׿�>doy?Nb�?�V�?��C?H��>���>��<Xh��-<=��>��>���>��;?l>�`?��<9w��KVe?հ?'�?>�u�>�κ>�_?�? �?E1?a0?>�?��?gg?-�?>y}?Lߒ?;{?E �?t#L?o\�?Z��?�w??�?L� >���?8g?D*?Z�?��>�� >�*>�o�>���>%8b='#�>za�>��e>`Є=��>�t>�Э>�o�>�P�>q8=2r�����<����Q[�Ea�^5�> �=�������yĿ5(�I��>��>鷋>³�>�(�>���=�_=�E�>���>�r>���?C`'?>��?�*?�>éD>��X?���?�W�?{>��q>���>�D^>��?��?"�<?(��?B ?P�?5>�?�s?(�}?0O"?���?�>�
�?�}?;�>�~�>��?C�>�/�>J��=ƞ'>D��>�ݦ>��>M��>E8i=�9���վ�e�������<P��J�y��:�hR�>��q%���L��b��씁������Ž�f:���>>�V[�'1���*�>�� >YJ�>�ł>�??�?�/>�G�>��F>�i�>�$k?	�4?�B?B>�݉>��>��>�	>�j�>�(k?��?���>h>|@�>�x�>Kk�>��_?�-?��?0#�?Mu�?F&�?O��?O�?-|[?TE/?4#�>�
�>�3l>��?�+?��?��>�{�>��>e��=印=N��=�b�>Q��>��>�)<WҴ�����a����)�exc�<e��0����=��+=�m�<��M=�������<�a=�hR=��>1�>0b�=(՚=�i�>�,�>�K+>��s>���?6�>��>?:��?F4�?$��?)��?.�?#�>���>�	>�<U7A=A��>�?4r�?�W	?���?��?�U>�4�>��;>�y>��t?�D?4�?jl?��? Gu?�?m�?8)?Tx�?I��?d��?v@�?cǩ?H��?:w@?V�?M"\>�C*=v��=ɹ>v(>07>���>�}�>E3U=���Ơ<�B��yNo=�'�=�*/=�ٱ>?��>��U>��>��=>�lu>%Ns=��:>��>�(>���>�I�>�L�>V��>F�>���?�8>��? �?tb�?�yI?��P?li�?*,_?Ty?{?�?>�a�>ò�>�$�>��E?�'8?�>�?�f?6�? yR>�zz>���>�8�>��>� ]>��y??2��?��>�II?t?��>��?�-?C^c?L�?Q؜?<.�>�,�?��?A�q?9��?��>_�;=���=R}�=Z_<�g��X���m��J�p ��ܤ�.Ӂ��	�=8�;�&ػ�1�=�6Z>!q�>�>3�>&ID=���>d�N>�u>�И?�G?#�?=�?-��?7�?$��?7�M?cN/?x?�||?�1?��V?��@?��>�қ>�כ?&K�??M>�tA>�n
���i�I��?���?jrA>�>J�b>��>>�0>�!|?T�>���>�-�?*�	?9`?�C>�3>�b>��?(=1?6�>��?6<�?U܋?�5?�?Sg�?�4>o
	>��=팦=�X�>E">O��>TIG>g�> F=i�=�<>6�>, ���׾#Y���>�����*���d>P(�>�~T>���>��]>�M�>�*�?�?��?��?:�?1e�?O?'7�?G,�?�3�?�?�cI?�([?2#?|�"?a]x?W��?�#�?a0U?7���I��?p��?��"?���?Q�@>��J?��L?��O>�2f?	��?;��?5�%?S�>�>�!�?�`?M��?673>��a>�QC?`L�?ZƠ>��>���?(y?GC?u*Z?yi?@xo?�|?B7>�֛>�/>���>E�,=�(=�(={w׾C��5��=�K�=�<��}�#��q>��p|�+����2�>�N�?Gޣ?V��?jW�?"��>�dI>�	�>���>x�L>��w?,w�?"\??��?C��?f�N?e��?xE�?�v�?�X�?���?����I��@#�?r�_?p�t??�?��>d��>�0�>�|�?Z9?H��?0�?^��?�v?J�i?-C�?(�!>��6?23?1��?+�?i��?8��?	H�?X�?o�?e9�?q?�?&ܽ?[{b?�6%?c�e?CF?H�k?G�%?:��? �6?, �?"*t>�ww>��=����3E��V�����W�> ��-^���>=qD?��>���?:�?f��?/l�>��=?	ۂ?��?)*?;�Y?@�7?G�o?N�b?<&c?@o�?N�?dw?�$j?��0?�3L?�]@�?l��?��?���?5ۼ?��?�A?'d?(a�?I
O?���?l F?Ȼ�?ԨI?���?*J=�p>3�<>M��?4�o?�Y?o�?��?�N?x�R?��Y?�P�?\_S?w�Q?b�o?^X�?��9?hg�?:�?Q �?d~L?e�?�DJ?l�?W��?`�?J��?6?]t>N��ha�K?<�vQ<��h���?"_?3C7?
d>�\�?(Ԋ??�C?��>�N>���>�?�;?P�?���?�n?ko0?b��?J|�? �?J^�?���?��?�S�?j�?�?��w?h�?V.�?��>/�>3�>�#>�d�?��?�ez?��?�ɉ?���?K.�?@�?v0�?hB?�? ?�Rg?���?�$�@ Tn?�?�?�qy?�w�?S��?5�\?C�o?8��?�B?cw>���;�j5��>�c�?1f?-��?1nS?Pb;?h]/?@.Y?��>�>�ٲ>�=>���>���?��?0�?� >��Q>�5>�q?	d?��>��>�$>���>죭?D(d?�?/�?Ry?��:?B"Y?
�	?2ɫ?��?�n?�WV?��?���?�DL?���?^��?.�b>�!B������>��'?+�D?~��?���?�U�?̑J?���?��-?�X�?�+?m&?�?��?�"�?ܥC?��:?H�?9�"?E?=_�?,/m=�t�:	�����پ`a����Կ2�$�����l�?n,?EB�?_=��Z��"�A��>��[?O�?ju?!�L?)	T?51�?9n�?94P?55�?��>�$�?
�Z?�	?&��?%��>�m>�e�?)'?Y%>.�<>��?0k>�P�>�Es>��?T�O?��D?�bf?��f?�~�?�0�?��1?�U?*�&>޹�>��=Љ�>ggZ?�?@�a?[�Q?��?��M?�oI?o��?s�n?���?�O@�?�H?�d\?|�?^�A?f��?F��?Z4=?Ltp?.�{?;?�?�w>��<'�=��j=�p�=m�*=�Y� &�!�־����>���?�|V?���?�Z�?seS?9��?&y?5�Q?R��?$��?�c?7�?{_?")0?3H�?L8? %�??�p?|�>ϰb>:L>�m�?_�>�VZ?+n6?A�?���?�(?�eV?���?��i?{-e?CU!>�x�>�e?	��>�>�b�?٦?$��?>i?\e�?��?��??�n?�)�?�}�@
p@$��@	��?�>?ǸG?|t�?B��?Pc�?] t?>E�?0�?��?֒?
=>�r>B=�{�=i �<2������&�";Y �=���>r�W>�ԡ?u�.?}��?%��?'�?/6�?"��?�?ڨ?'u??)��?;�?*�?�:?�K@��>��h>�W@>�f�>U�<�6S�6Ŀ�<����>�S[?^�?rvL?|�a?q}�?n�?�<`?Uw=?	��>��>�ϩ>�~V? �>�2>�?
�h?)#�?I�e?H!�?Yɮ?�,S?���?�K?�L�?Ѭ�?ڍ�?�X�?�se?��>?X�?\�s?=&�?(=?7�L?��>���>��>�z>�b�> ��=p�=���;Wr>�?=��R>%rn>�ߺ>��i=+�=�1>��$?"v�?6��?2x�?(rY?(8�?(�4?PU? "�?,�о֏?�0>�;6>�U&>s>T�7>x��M�p%>�.��q7/�cM��ǚC����>x��?>��?`�??	?�?
��?q� ?k޶>�3b>��q>פ(>�v >��D>x�>h�->��.>�"�>��{?$u�?k-�?�T]?��t?�'?�N?_�?:�?
>��?}?6Z�?i��?7��?8�?:�?Fح?��>�_=�n�.��h�T��Y��>J*>ȱ]>��|?),�?(�z?E��?|��?VT�?0��?@�`?F[$?Pf�?Rb!?S1�>�Α?9u�?��e?%r�?�>ڻ>��>@>nֲ=�Y��=^�>ղ�M���9���y3�i�K>�,{>�d?�%��@=�	�=c��>���>�l>\�U<�L)=�/=�PG<E��<{8�>6��>�s�>��8>��+?%�#?g�q?��?��?yњ?I��?:�??b?>�l�?�?IZ?+2?G��?�uB?g�F?>%�?ق>�FY;��L=9��=�MC=_:D>��>��^?;�o?d��?/i>�Z?=b�?[��?U=6?#`?*�p?0�l?�?4y>T��>`�E?7��?��>١�>��L>�p�>��a>A{=������I��>��n>Y/���,=\��>�k�?M!=>�y=!��>]�@=��2>cb
>\$Z=��7<2��yF�9¾!:�1��ǖH���n;��>q��>��U>��?�8?3�w>��>���>�4>��>�Xf>M��>F$>X�C>�+�>��>�}'>�'�>bI�>碮?zUT?gK�? �v?��?D�?���t��?ǳ?< ?;��?�?�$�?�u�>�i�?/T�?66�?W��?�}?#%�? 8E=�H�?�"?Z�>�`�>�sK?2m�>ɍ�>���>/�=;�>���I���I���I���I��?�k,?�F��?r�?T�V>j k>�_�>��>MS=}i����u/&����Q���M�k<�`�=��>=��>Tͥ>���>[�>��<>��>��?-8>��H>�*m>�ү>��,��f6�t�����>��>�E��I��?ɼK��I�������/�����'��2�*�n���$�>}����:w�W+�=���I���I���I��?��?����-[�?f��?�%���꾨�0�hu>�{>�0�>���>�=>���?}z��I���I���I���I���I���I��?2F�?�����޽�l�>��r�*ɯ=�F��SZ���
>fV�>w����(O��EW��+0�kFN?�g?#>�5'?�R3?��>ѹ��>5��,���0M?�i@N�@�I���I���I���I���I���I���I���I���I����A������|�j�I��    �%Ǘ�I���I���I�ʷ������!�z�K?SG��I���I��>�#?g0)�I�ʾ��h�5��=.l�?f:�?S�]?_$(>�RV�I���I��>��>���I���I���I���I���I���I���I��?�r�?���?�c�>�n��I���I���I���I���I���'��I�ʼ�I��r��SS����J��m����W�)U8�I�ʾ��j�I���'[T�I���I���I���I���I���I���I���I���I�ʾ�u>+�>>���I���I�ʿ�����L�	N��I���I���I���I���_��e��_��_
��I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I���I��>�p#<5�s�.i��I���I���I���I���I���I���I���I��