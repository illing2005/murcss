CDF       
      lon    H   lat    $         CDI       GClimate Data Interface version 1.6.1 (http://code.zmaw.de/projects/cdi)    Conventions       CF-1.4     history      VWed Jul 09 12:03:24 2014: cdo sub /tmp/murcss/cache/constantField.nc /tmp/murcss/cache/2_5_tas_miklip_initialized_mpi-esm-lr_initialized_1960-2000rmss_sqrt.nc /tmp/murcss/output/20140709-120305_tas_miklip_initialized_mpi-m_mpi-esm-lr_initialized_miklip_uninitialized_mpi-m_mpi-esm-lr_uninitialized_HadCRUT_1960-2000/2-5/miklip_initialized_mpi-esm-lr_initialized_input1/msss/2_5_tas_miklip_initialized_mpi-esm-lr_initialized_1960-2000_rmss.nc
Wed Jul 09 12:03:24 2014: cdo sqrt -sub /tmp/murcss/cache/constantField.nc /tmp/murcss/output/20140709-120305_tas_miklip_initialized_mpi-m_mpi-esm-lr_initialized_miklip_uninitialized_mpi-m_mpi-esm-lr_uninitialized_HadCRUT_1960-2000/2-5/miklip_initialized_mpi-esm-lr_initialized_input1/msss/2_5_tas_miklip_initialized_mpi-esm-lr_initialized_1960-2000_msss.nc /tmp/murcss/cache/2_5_tas_miklip_initialized_mpi-esm-lr_initialized_1960-2000rmss_sqrt.nc
Wed Jul 09 12:03:23 2014: cdo sub -sqr /tmp/murcss/output/20140709-120305_tas_miklip_initialized_mpi-m_mpi-esm-lr_initialized_miklip_uninitialized_mpi-m_mpi-esm-lr_uninitialized_HadCRUT_1960-2000/2-5/miklip_initialized_mpi-esm-lr_initialized_input1/msss/2_5_tas_miklip_initialized_mpi-esm-lr_initialized_1960-2000_correlation.nc -sqr /tmp/murcss/output/20140709-120305_tas_miklip_initialized_mpi-m_mpi-esm-lr_initialized_miklip_uninitialized_mpi-m_mpi-esm-lr_uninitialized_HadCRUT_1960-2000/2-5/miklip_initialized_mpi-esm-lr_initialized_input1/msss/2_5_tas_miklip_initialized_mpi-esm-lr_initialized_1960-2000_conditional_bias.nc /tmp/murcss/output/20140709-120305_tas_miklip_initialized_mpi-m_mpi-esm-lr_initialized_miklip_uninitialized_mpi-m_mpi-esm-lr_uninitialized_HadCRUT_1960-2000/2-5/miklip_initialized_mpi-esm-lr_initialized_input1/msss/2_5_tas_miklip_initialized_mpi-esm-lr_initialized_1960-2000_msss.nc
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
Updated at 26/06/2014 15:02:28     CDO       HClimate Data Operators version 1.6.1 (http://code.zmaw.de/projects/cdo)          lon                 standard_name         	longitude      	long_name         	longitude      units         degrees_east   axis      X        @      lat                standard_name         latitude   	long_name         latitude   units         degrees_north      axis      Y           `   const                      
_FillValue        ����   missing_value         ����     (�  �        @      @$      @.      @4      @9      @>      @A�     @D      @F�     @I      @K�     @N      @P@     @Q�     @R�     @T      @U@     @V�     @W�     @Y      @Z@     @[�     @\�     @^      @_@     @`@     @`�     @a�     @b      @b�     @c`     @d      @d�     @e@     @e�     @f�     @g      @g�     @h`     @i      @i�     @j@     @j�     @k�     @l      @l�     @m`     @n      @n�     @o@     @o�     @p@     @p�     @p�     @q0     @q�     @q�     @r      @rp     @r�     @s     @s`     @s�     @t      @tP     @t�     @t�     @u@     @u�     @u�     @v0     �U�     �T�     �S`     �R      �P�     �O@     �L�     �J@     �G�     �E@     �B�     �@@     �;�     �6�     �1�     �)      �      �      @      @      @)      @1�     @6�     @;�     @@@     @B�     @E@     @G�     @J@     @L�     @O@     @P�     @R      @S`     @T�     @U�     ��h ��@����������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������>2w�>4*��������������������������������������ρ̾g�}�@����j�J�cEԽv��<mI�;
T �t�N��̾2]8����������������=XYp��� ����������������������������������������=/� ;݄ �zJ��i�$�0 ����Ƚ��`����>^���������������������>��>�P����������������������������������������������������������������������������������������������������������������>W)�� ��נ�;V��_$��������\�`��=Q�`>x<��/�P����������<����������俩�,�^&���b:�}�h�?}t����'�������������>�xh�����wmh�Y�[6t��=t�v�侟���;��I頾�x����=����o��`�> ��>�+�>���<lv�<
 ����U@��>Y�@>q{��%v�>�C�>]bh��E�>?�>W(��u�>}�>Qf@>Z��=��8=�c8<�� �VD>|���V ���U�ȿ�wX�Q,.� :����0��¿i�V�!��=���?��>�;l>"z�� ��>��>!6�=Hwн@������������:��˒��e��[�ؿT������*��<���==D<�J`��P>��D>��>�袾6FP���*��0��+���8� `��p�7�H�k¿9솾�о1���F�����ؾ��E �K:@��7@�n� �6��lo`�����ػ�a =��=��X�k�@� (�<���=Cd�g����Q�>5X`��F >7b�bB���- ��%��A6@�������XO`�C0 �Ҡ��O =���<��=֧@=��(�!z���(��y@=U`<���̀�.[`�&9 �����������l�y�>/�����<�"�=�0>k��>�*>�k�?�x>��*��f���{ ��V<���t���\��,Ⱦ�"��ǰ���\�`&���0�$i��x�r� ��:��9��݀�z�Q����P�,�x�S�P�AC(�� >!� =i-0>.J�>$5��[�P=�0��‼9|�=ی=�_h=���<ҋ���� =׆P=��>RT=�E0=�I`��)л3� ��v��X���,����y*x�'���������������
GP�=о�ؽ���<+����;�Z�X��� =��=�i�>V��>���>�9�>BM$;|� �ΰнݾp�n�H�-$��%j�[���>�V�@�{ �]�h�\5ؾg>������Df��yp���=?$0�����0���о
� ��Ɣ���|H�ͱp��ɀ������@=X�`<V �􇀾v�h����=�S >%l>,ܐ>+Ո;�} ���=��=9�=�z@=���  �UW�.�(��� ����=��`�9BX����ǻ���� �����'�x�p�����0���P�S�`������S`�y����?�=�nh=�{�==�b(>S��=��>$Ը>��>?��;1� <]����Mо��X�:(�0��:� ���(�q�X����"j�����J��o�='p �] =T0�>��@>�U�>�n�>�T�>�"t>	¤=w����ཧ���:<�=���>:>�̾3��q%`�UHؽ�t��y@;���1��;Z ���侓<���0�Q� =~��=����p�9\�^q���(�=����;���� <
��=�\X=&�0��} �B���A����T�2߈��|��IA@��� �u� =��0=_[�=4 =+��=o׀>	���聴>2�\=�a�<�� ==��<`P�=ؓ0<���3 ��`���P��������H�8�Zw`������>W�L=���> �\><��>1�t>W�>6�>�&>צ&>��>��P=�j�>u��>̎$>`&�<�`@;�� =�i8=�U���<@>�Q>:`<!L@�����ecH=�Ĉ=ޏH=�Wx��؀�n%���S����ԾSp��(�>>(U�=�g �ᚐ��򠾃��uU`��^�q���p��DF��E ؾȽ��0=��>K�\>%��=�� ��I�<� �=�޿z:b=�
 =���>M,�=(�p��G@=�W��%yཅu��@�%����2t=�
�<-y�>���>��b>ᬼ����m��=�٠=��>Y�t>��>�\�>�`@>���>��^>��(>v:,>vf\>��>��F>�U`>�'f�Ў ��꾻E(��ؾ�簿$w��F�ξ���>��H>��=���<�@@<2R �.�
�༼����1 <�f�=���>�>W=��=�x`=|s ��@�+��=Y������5�0���ľ `�9����7��+(R�.^�7�c�@=o���Ȋ=�E`>E�l�Y�������Qн���������85 ���̾;?h>���?*�?�J>���>���>J��>� ?ݺ>�#>أ�>�j>\2L>�=
�@�Ud�C >&,>��>�76>�p$>uX�>� ���ؿ	T�Ѹ(��쿠�Կ�!辑�=I�`���P=�͘>7`=���=�r(�0&����p;n >�>p�l>L��>"�> i�=R��D�@�H�;����,e`��f$�9"ྊ�о��(�YFl��%����o��<�@�ux����+R�����>��>>�P���P�Q�@��� <�W ;� =�8�=���>[>��6>gIH>��?��>@x�� �a <T�?�*?!�t>��N>�|6>;�=��x��h ����>xt>�Z?ܕ?�@>�"��gϸ�;s�����;J��% �%�v�⠽�s0<�X@�>ܰ��������
0>H��>/�h=$��>8�T>�R.>=\>n@=��>.L=�x�=e�=�0=f��=;  ��YP�;� >e,�_H��j�`�UN��4�1���\�U������#�D��!�� �!�z>�\>��6=G`�>&��>�TN>�??�i?��? -�>��v;#� =?-0�x� >�t�"AཥlP�n* >��(>�G�>���>p&�>���>��j>���>���>���T0X�*4Ⱦ����,��`ྵ�̾����4X2���ܾ�8� 'ʿ>��.w���	����`��j�>L=�>U�x>4�>��$>=��><�=g� =<�>4=�Wp=��><���` >��>�k�>�w�8 ��6 ������R<�N��U檿b����(��N<��4���R��� ��y�>�X>��>���>�b>��J>m�=�y(��Xp�#v��� =�J �8�ȿ����^>�E�>�Dh>�|�O�`>T����>;�>f��>��>���>`?��(����:�b��|��(���:�=�`�T�����2�g�t�z���C��`о�3�� >w�>�zĽl@�����>#�>cB<=�t�/ <��@=8�0>N|>�O> ��[�����7���� ������D�^��}V�K�n���Qل��Tľ0h�[���D꾿+伃~�?	2^>��>�����r��?J��� ��K��G侞�=�p��͜�l
 ���%��?8�>���>��>QO�>�Ҿ���A�Ƚ��о"Ř=��>	�<Ui@>�L>P�>Cո=笨<�- >�m�>�Ø>l�8��O���帾�� >�>�Ⱦ�s<��@>���>��0�oW@��TP�/� �pI � �>���\ĽS=`�* 辢+併�@=�-�=���������P���<�������辛��4��f`������V���F`�����0T������g>c�>� >��>��n���ؿ2���B�@�Y:��Z�R������,4�2�`=��X��İ=��8��8��)�V=�
0>J������W��3 =�̀;�>\I >7}�>Gp|>�X>7�>P��>���?(u>��j? C?:��? ?�?��?i�? Ǥ>��"�S8 >��=��>&d�=��н5� �u� �D���xy8�ǂ4��?�Me`�ٷp��о�騾����`�U%�_B ���,��4������\��������i*����Q* ��Wh������8��t�����	� ��`�;<(�ྵ-��a�:��п5�0��9���?�>�Q�>�
>�=�b��Ԓ���>���>�"�>���^��������|=�fp>�k>��>�T@>�>��=��@>8�>�H|?>>�>x#���>��<>�AƼ�������=��qr��k���n��{t�Lr��-���e���Jh��Z������6$��@���b�F�X� �к�)���#TP�RB���-`���t���D���T�ɣP�!�@��΄��8�Ӧ��>�ο�z���>�=��P���V{���+��_H=�� <���=��о�T8�)�p��+�!��
P���? o�>�Xf�%�ƿ)Y2���v��H���=Am`=�<P>�>��f>�"�>�H>��>��?x>�j�]��fg��Pp>p� >	T =�1���C�������a��L8������0��A���	u��2�h���|��c��߈̾�Gd����_��Vl����kT�|N�ȭ��_�̿�0
���J�ED|�����ɴ���l��^(�c+��},`�^���%�����>����ఽ�9�>���>��>��d>��">a#�<?��
I����о�	�+������\ڐ>���>�J�>���>[�l=� ��1��'&�$�V�!��>ckT>��>>%1>d�;$� =HX�>X� >�^�>�M�>��j? �b?e2?(��?�?�?T? ��:S��K������x��Q���I�������qP��}�Ҟ�� $��:Ծ��x���H��k���b��jBо��(��Rd�mdD���P�Zմ�� ����ľt\x�Ƣ@�aK��V�:8`=�h����>5��>⬈>渆>�x*?�>��>��(>z�P>�p�.� ����H���L>���>Ќ�>�q">��`;�	 �]������GX��Jp����9�(>-�>�Q�>8ٌ�K�=�k�>m���%ް=���>�. >�l>�Ub>�W���0>��>�ͼ>�4z<�Q��Bȿ�� �ֿ2�"�J�(���2�ڿ>!̿(YB�
N����4��z4���0��%��p�úh��пk�����(�����#�,���,���t=��>wK�=��>���>��X>4Mh>Xx4>���>���>���>���>��?5�=)��;4��7 >���>�`��u���x��� ����?	�?)
��CJ ��}��X����~�<�/གྷ������>R� >à*>4�ؾ�lt�7",�g�>�
>o3�<6 >n4>�^>���>n�?-�>(�Ѥ��P�M��R��xR�ƆԾ�_4��@�� �=���1(��?�M�لԿ
�¿�Z�I<��MH����̢����_��8*���h<���0<W0�>VL<�� >1W(>y�4>�t>'�0>��(>�h2>���>��>�>�in?+��>�3�? s?��>꬚>\�����>x�>�Y�>�A">��0���>�%�>�����=�R>�i�>���=[��������8=�>�^�>{����T��>�?I>�i��b> �U�>��>�[�?�?9�>�Ͷ=�� <e�@�3��ݽ4�Q��>~��Rޢ�+� �f�I<�B�J��el��	̾�kX� �<��h��Ծ{�ؼ� �>��V?_�>Ǜ>Z�D�*%`�%4������9��M��>��Z>3�$=�?P>��l>�i�>���>��F>Ƽ�>�Q?/��?oH����>�~�>�kN>�h"=�_�>�<���(������� >�;�>�#�=�R�>��b>�$>��X>�%>dE���0=��`>��>��?(P>��B=��p>�k�?��>Şf>9\>O� >��?�?��>��>��j>0>��f>[O�>{\>5�P>"<�A��:p��>������^���J��������Ծ�b;�	�����>�y`?'i�>�z`�g�P=c�`=�Ԙ>T�>��>�>�
\>�}�>�6>��>���>�,b>�ǆ?h>�g�?��>�@�>�Q�>˚�>피>���=�"8=��=qU�>��>���>�a�>��h>�~>�m�>�+�>�5~���J������j>���<���>&��>�(r>ݎ?�J?"�?�c>�p?t�>�k�>�>��>ȷ�>��>�`>̹�>��>��>�3�>z��>f�>��=��=���RX���`������j��� ����>,��>���=���T%�>���>�D>>A���T�ཱུ̇�P=�%�>N�D>ۛ�>�WJ>؂�>��V>�SX>��>��b>��p>�R>�4�>�B�>�J�>�ަ>��">�,<Ǵ������נ�W4 ��� =��`>�I�>��>�S.>��>���>�?>>Ʒ�>�d�?+8?-�">�wd>�
�>�F>؀>�,?��>�h0>�:8>U�H=��=��P=>=辒���mн���=���=Ɓ�=�x>h�>� =��(<;F���z ��.`��� �:�8��| >(�>��>_[� ������ɒ�=�7(>���b� ��8��@��CP>�}h>˝�=��h>@h>҉�>!+<=(�0>?� >�Q�?�f>�>�EH?PH?5�?��>��>�J��=``�x(`���(��P>yL>2d>؏,>�l>��@>ߩh>���>��>�kX>��=DA >�2�>��l>�|*>��j>��>`�D=��0>���>9L���ֿ��ȿc���T� �k[��6������.@=��=��(<���Ͷ �/@�KH�iπ=��@>�=�>!��>[G >��4>T��>~9`>O��<�����}P=w��=�O�>Bv,>T%(�g! ���>j��>���ȗp����<P	@�^ �H ���>��?O�>��>В�>���>��Z>�m>�02=𮠽�8@����&¾�� =��>ZS>I?�>��p>���>i�`>G8>�<>F��>�hB>�� >���>���>���>gh�>,Y>�x>Bc�>�ǰ>��z=ݲ�<�/ ==`���t�G���t����Ԋ̾�)��������l�ߪ(���,���>���>�V�>�`>Ʌz>�>:C4>c,�>h[�=ʳh<-S�=G��>��>4/�>wDx>
�@:� =��>2��=�Ĉ������Bؽ�� =xⰻ�N >��V>���>��>��>>���>V?�>L�4>L� =�J��=m ���@<�݀�ý@�h@=l�>
>_��>��>�]�>�]�>f�<>޴>�d>L�$>d�0>��>�s>�z >*�>��>L�>�-�>���=�?�=8p=1��=�� ���h��]��IcH�y�t�t�H�Xn��)
¾��X��1��<Ϡ��?�>M�(>���>!�>]�8>�G(>3{�=�� =��>�,>,D=��=!��;���<��`>���豈��H�@P����1@��`辵9ྐ,4���>(И>�(�>|4>M��>5~>W+ =�cp<vX@�7�འI ����=z@�h‽'[�=�0>�8>J�p>+q�>&<>3 >x0�>]7T><D>%��>"Dt><0H>L��>/�8>k`>10�>%��>�X>��=�nP�����M0�� �l^X������S��~���}�L�!�n��T �9z����� �a$о�Od�&�X�9��>^��>��v>�9d>Y�0>/s4>	�=�K@=�.�=��轩 =�����( ���p��~���B�����M�������=D�Ƒ`�q���`�3���w� =Ò�>/��<ez�<�K�=���=�@���@�xK �>�`��g0�+c�JN��\Fh� ���@��� =�@>�`>McL>�i�>|�0>>�l>J�=�lx<ڲ�Yt <�� >�wr>��>�J�>#�>�FN>�u�=۸h�	�R�K �12�E����7��"J�d䠽�  �)�=���=���>R0�>�)�>��@>G^\>���>�t�>�$>P� >�0�ӊ�>6A�>�d>� =�����l0���$��=��X�*��]��m���\�����1<�� ��n�:�(�\��nY@<O�@�Fb���H��$�2Ƚ�&@�"j8��L��Ϥ��C�������p��x]���倽W =���>V�>�Q�>��>�-�>IB�>>��;f^ <t?@�u��=�`>���>-�>xd�>��h>���>?�=Dཡ�@�;B��ؾ&����I`���L�
S >�>��>b$t��� >n�\>��">�w�>[�><�P>JA�=�҈<�H �Ga��<��>��l>+0��о��P��8T�t~��Ƨܾ񩴾��$������G@�UO��e� ��] ����>�^J��t��)3��y�п�B��
�r���=Kƿb��n���9̿������l�1�|��ʾ� �M<�����c��=_>�>(�����@��� ����̸�x����о��,���(�mpP��5���<����P� ��(� =�e�=�PX:J  <���=&� =
9���;�m >��j>� L=�>�t*>�zD��`>NP>gJ�>�j�>���=��p:Z � a=>)��>	�����t�&�>�^8��a �K�EBX�\[ؿ�T����������������>ܟ�>O�����>0�>�uV�w�N��<��r��-��J�ʿo3r��*��WZ��FL��������v�������`�p���h���p�8�@��6�>Yt�*�`�>�཯H0�(�P�ҕ4��G����\��<`��A ����;Z� ��� �����	>��p�6�,���4��Q����.0�� ��������I��)[�������������>��>������=���>;�8��Z���о�0 �����J4��h��- �
3X=�q0������������������������>���<s���z����`���� ���h�N�Z�G)8� ��J俩N���{B��������A�=��=D��C�H>��>�V����P���8���P=��>����������������������������������������t�'� ���p���������Ũ��������������
9�Mz`��b�����?IЄ�����������2>������O�<��g��Ţa>�5�>$�X>/���E�����������s � �����������������������������>�pX>�O�>�p���v�����������������������"������.��ی����B��R&p�1����� �����>��������������������������������������������������}�d��*���������li��}����������������������]`���辒�侒� ��������������������������������������������������������������������������˜°c���Z��������������������������������