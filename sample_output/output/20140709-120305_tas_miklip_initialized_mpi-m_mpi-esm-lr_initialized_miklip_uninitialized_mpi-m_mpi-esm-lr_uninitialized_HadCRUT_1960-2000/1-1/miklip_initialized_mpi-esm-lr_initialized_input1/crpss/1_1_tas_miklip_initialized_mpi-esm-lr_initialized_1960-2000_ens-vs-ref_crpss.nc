CDF       
      lon    H   lat    $         CDI       GClimate Data Interface version 1.6.1 (http://code.zmaw.de/projects/cdi)    Conventions       CF-1.4     history       �Wed Jul 09 12:04:32 2014: cdo sub /tmp/murcss/cache/constantField.nc /tmp/murcss/cache/crpsSum11_1_tas_miklip_initialized_mpi-esm-lr_initialized_1960-2000_ens-vs-refdiv /tmp/murcss/output/20140709-120305_tas_miklip_initialized_mpi-m_mpi-esm-lr_initialized_miklip_uninitialized_mpi-m_mpi-esm-lr_uninitialized_HadCRUT_1960-2000/1-1/miklip_initialized_mpi-esm-lr_initialized_input1/crpss/1_1_tas_miklip_initialized_mpi-esm-lr_initialized_1960-2000_ens-vs-ref_crpss.nc
Wed Jul 09 12:04:32 2014: cdo div /tmp/murcss/cache/crpsSum11_1_tas_miklip_initialized_mpi-esm-lr_initialized_1960-2000_ens-vs-ref /tmp/murcss/cache/crpsSum21_1_tas_miklip_initialized_mpi-esm-lr_initialized_1960-2000_ens-vs-ref /tmp/murcss/cache/crpsSum11_1_tas_miklip_initialized_mpi-esm-lr_initialized_1960-2000_ens-vs-refdiv
Wed Jul 09 12:04:31 2014: cdo enssum /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1985_r1i1p1_198601-199012.ncGPOU5Q.remapped_selYear_1-1293973214METimmean_ANOMALIE_BIASCORRECTED2_ENSMEANcrpsref1985 /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1990_r1i1p1_199101-199512.nc4K2BBX.remapped_selYear_1-129438VA76HLTimmean_ANOMALIE_BIASCORRECTED2_ENSMEANcrpsref1990 /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1960_r1i1p1_196101-196512.ncIUUWGS.remapped_selYear_1-1291967HA85UTimmean_ANOMALIE_BIASCORRECTED2_ENSMEANcrpsref1960 /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1995_r1i1p1_199601-200012.ncQ92ZO9.remapped_selYear_1-129478S5PPGZTimmean_ANOMALIE_BIASCORRECTED2_ENSMEANcrpsref1995 /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1965_r1i1p1_196601-197012.ncVHO52I.remapped_selYear_1-129236C1CNJ9Timmean_ANOMALIE_BIASCORRECTED2_ENSMEANcrpsref1965 /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized2000_r1i1p1_200101-200512.nc9ZZXDF.remapped_selYear_1-129519FMKE6LTimmean_ANOMALIE_BIASCORRECTED2_ENSMEANcrpsref2000
Wed Jul 09 12:04:23 2014: cdo ensmean /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1980_r1i1p1_198101-198512.nc9CAON1.remapped_selYear_1-129357P7C9SGTimmean_ANOMALIE_BIASCORRECTED2 /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1980_r2i1p1_198101-198512.ncWVL82Q.remapped_selYear_1-1293583MLOU5Timmean_ANOMALIE_BIASCORRECTED2 /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1980_r1i1p1_198101-198512.nc9CAON1.remapped_selYear_1-129357P7C9SGTimmean_ANOMALIE_BIASCORRECTED2_ENSMEAN
Wed Jul 09 12:04:23 2014: cdo add /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1980_r2i1p1_198101-198512.ncWVL82Q.remapped_selYear_1-1293583MLOU5Timmean_ANOMALIE /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1980_r1i1p1_198101-198512.nc9CAON1.remapped_selYear_1-129357P7C9SGTimmean_ENSMEAN_ANOMALIE_shift /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1980_r2i1p1_198101-198512.ncWVL82Q.remapped_selYear_1-1293583MLOU5Timmean_ANOMALIE_BIASCORRECTED2
Wed Jul 09 12:04:23 2014: cdo mul /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1980_r1i1p1_198101-198512.nc9CAON1.remapped_selYear_1-129357P7C9SGTimmean_ENSMEAN_ANOMALIE /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1980_r1i1p1_198101-198512.nc9CAON1.remapped_selYear_1-129357P7C9SGTimmean_ENSMEAN_ANOMALIE_bias_shift /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1980_r1i1p1_198101-198512.nc9CAON1.remapped_selYear_1-129357P7C9SGTimmean_ENSMEAN_ANOMALIE_shift
Wed Jul 09 12:04:23 2014: cdo subc,1 /tmp/murcss/cache/conditional_biasNew.nc /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1980_r1i1p1_198101-198512.nc9CAON1.remapped_selYear_1-129357P7C9SGTimmean_ENSMEAN_ANOMALIE_bias_shift
Wed Jul 09 12:04:22 2014: cdo mul /tmp/murcss/cache/_correlation.nc /tmp/cdoPyXGjORa /tmp/murcss/cache/conditional_biasNew.nc
Wed Jul 09 12:04:22 2014: cdo div /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1961-19653G369J.remapped_selYear_1-126541JAKOW5Timmean_ANOMALIEmeanAgain11_variance /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1960_r1i1p1_196101-196512.ncIUUWGS.remapped_selYear_1-1291967HA85UTimmean_ENSMEAN_ANOMALIE_variance /tmp/cdoPyXGjORa
Wed Jul 09 12:04:22 2014: cdo timstd /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1960_r1i1p1_196101-196512.ncIUUWGS.remapped_selYear_1-1291967HA85UTimmean_ENSMEAN_ANOMALIE_merged.nc /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1960_r1i1p1_196101-196512.ncIUUWGS.remapped_selYear_1-1291967HA85UTimmean_ENSMEAN_ANOMALIE_variance
Wed Jul 09 12:04:22 2014: cdo mergetime /tmp/cdoPya4zqTz /tmp/cdoPyC1e2_z /tmp/cdoPy5tsg5L /tmp/cdoPyzJIYcB /tmp/cdoPyotfBA9 /tmp/cdoPyLuGzAE /tmp/cdoPyRD_pLU /tmp/cdoPyjADuis /tmp/cdoPy2TMetZ /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1960_r1i1p1_196101-196512.ncIUUWGS.remapped_selYear_1-1291967HA85UTimmean_ENSMEAN_ANOMALIE_merged.nc
Wed Jul 09 12:04:21 2014: cdo yearmean /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized2000_r1i1p1_200101-200512.nc9ZZXDF.remapped_selYear_1-129519FMKE6LTimmean_ENSMEAN_ANOMALIE /tmp/cdoPy2TMetZ
Wed Jul 09 12:04:20 2014: cdo sub /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized2000_r1i1p1_200101-200512.nc9ZZXDF.remapped_selYear_1-129519FMKE6LTimmean_ENSMEAN /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1985_r1i1p1_198601-199012.ncGPOU5Q.remapped_selYear_1-1293973214METimmean_ENSMEAN_2000_model1-1FD0E4Y /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized2000_r1i1p1_200101-200512.nc9ZZXDF.remapped_selYear_1-129519FMKE6LTimmean_ENSMEAN_ANOMALIE
Wed Jul 09 12:04:18 2014: cdo ensmean /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1985_r1i1p1_198601-199012.ncGPOU5Q.remapped_selYear_1-1293973214METimmean_ENSMEAN /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1990_r1i1p1_199101-199512.nc4K2BBX.remapped_selYear_1-129438VA76HLTimmean_ENSMEAN /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1960_r1i1p1_196101-196512.ncIUUWGS.remapped_selYear_1-1291967HA85UTimmean_ENSMEAN /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1995_r1i1p1_199601-200012.ncQ92ZO9.remapped_selYear_1-129478S5PPGZTimmean_ENSMEAN /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1965_r1i1p1_196601-197012.ncVHO52I.remapped_selYear_1-129236C1CNJ9Timmean_ENSMEAN /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1970_r1i1p1_197101-197512.nc8W2Y2K.remapped_selYear_1-129276WQJFFLTimmean_ENSMEAN /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1975_r1i1p1_197601-198012.ncXW42Y4.remapped_selYear_1-129316ICW30RTimmean_ENSMEAN
Wed Jul 09 12:04:17 2014: cdo ensmean /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1980_r1i1p1_198101-198512.nc9CAON1.remapped_selYear_1-129357P7C9SGTimmean /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1980_r2i1p1_198101-198512.ncWVL82Q.remapped_selYear_1-1293583MLOU5Timmean /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1980_r1i1p1_198101-198512.nc9CAON1.remapped_selYear_1-129357P7C9SGTimmean_ENSMEAN
Wed Jul 09 12:04:14 2014: cdo timmean /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1980_r2i1p1_198101-198512.ncWVL82Q.remapped_selYear_1-1293583MLOU5 /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1980_r2i1p1_198101-198512.ncWVL82Q.remapped_selYear_1-1293583MLOU5Timmean
Wed Jul 09 12:04:14 2014: cdo seltimestep,1 /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1980_r2i1p1_198101-198512.ncWVL82Q.remapped29358yearmean /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1980_r2i1p1_198101-198512.ncWVL82Q.remapped_selYear_1-1293583MLOU5
Wed Jul 09 12:04:14 2014: cdo yearmean /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1980_r2i1p1_198101-198512.ncWVL82Q.remapped /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1980_r2i1p1_198101-198512.ncWVL82Q.remapped29358yearmean
Wed Jul 09 12:04:05 2014: cdo remapcon,r72x36 /home/illing/workspace/murcss/sample_data/miklip/initialized/mpi-m/mpi-esm-lr/initialized1980/mon/atmos/tas/r2i1p1/tas_Amon_mpi-esm-lr_initialized1980_r2i1p1_198101-198512.nc /tmp/murcss/cache/tas_Amon_mpi-esm-lr_initialized1980_r2i1p1_198101-198512.ncWVL82Q.remapped
Tue Jul 08 16:59:36 2014: cdo remapcon,r72x36 -selyear,1981,1982,1983,1984,1985 /miklip/integration/data4miklip/model/global/miklip/baseline1/output/MPI-M/MPI-ESM-LR/decs4e1980/mon/atmos/tas/r2i1p1/tas_Amon_MPI-ESM-LR_decs4e1980_r2i1p1_198101-199012.nc /scratch/b324031/CMORSTRUCTURE3/miklip/initialized/mpi-m/mpi-esm-lr/initialized1980/mon/atmos/tas/r2i1p1/tas_Amon_mpi-esm-lr_initialized1980_r2i1p1_198101-198512.nc
Model raw output postprocessing with modelling environment (IMDI) at DKRZ: URL: cosmos/tags/mpiesm-1.0.02, REV: 5969 2013-03-18T11:06:04Z CMOR rewrote data to comply with CF standards and MiKlip requirements.   CDO       HClimate Data Operators version 1.6.1 (http://code.zmaw.de/projects/cdo)          lon                 standard_name         	longitude      	long_name         	longitude      units         degrees_east   axis      X        @  #d   lat                standard_name         latitude   	long_name         latitude   units         degrees_north      axis      Y           %�   const                      
_FillValue        ����   missing_value         ����     (�  &�        @      @$      @.      @4      @9      @>      @A�     @D      @F�     @I      @K�     @N      @P@     @Q�     @R�     @T      @U@     @V�     @W�     @Y      @Z@     @[�     @\�     @^      @_@     @`@     @`�     @a�     @b      @b�     @c`     @d      @d�     @e@     @e�     @f�     @g      @g�     @h`     @i      @i�     @j@     @j�     @k�     @l      @l�     @m`     @n      @n�     @o@     @o�     @p@     @p�     @p�     @q0     @q�     @q�     @r      @rp     @r�     @s     @s`     @s�     @t      @tP     @t�     @t�     @u@     @u�     @u�     @v0     �U�     �T�     �S`     �R      �P�     �O@     �L�     �J@     �G�     �E@     �B�     �@@     �;�     �6�     �1�     �)      �      �      @      @      @)      @1�     @6�     @;�     @@@     @B�     @E@     @G�     @J@     @L�     @O@     @P�     @R      @S`     @T�     @U�     �yv�u�������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������°�� �������������������������������������Xx :�� <AL ��[@�H#p�_@�W@�SJ�1�<� =ޣx����������������������� ��� ����������������������������������������������= <�a�f� �f� ��@�ޜ �>� ���P��������������L@�����3�8!@ ���������������������������������������������������������������������������������������������������������������������r ��0 �� �Q� <t"����@�	�u��[���4��z���7 ���&jX�.��� =���=��=��<%,�<����.�=i0=;�;s\ ������������������������=}�@�������t��.t���`��}���r �k� <A,����@��@ ��3 �k� �� �����N���н�6��)K���@�Ӑ�ˏ0�Zw�ʩ�/g ��J@��s���@ �=�`;� ;�"��*� ��P�]_ ��� �-. ��r��=༩:��ʓ��h ��W ;A; �����J� ����� ��A �� ��� ���𽱖`��`���о@y�J8�U�:H� �MH���� ;� �����Q���O@��< ��нֳ ��� < ��� � �� �Z  ����;�U ;�E������r� �!�������q �������@�'7��O� ��� �ߐ �I���g@86@ �[, ;���~� ��� :�z ;����	�����лL ��Y �qX �xx`�Y� �Ӿ���`;� ��H�:�l ��� �4��ޠ�����Ԇ@��[�:>L ����x��7T�P������^��j��B`���x�ʽ`�̅�Ep�^�H����F� �a�辄$����c���$���P4о�Z`���0�G@�O� ��v �����ō�����<#[���� ;������ ���`��� ��� ����&}��� �&��;h� ����U���Y�� ��l`���� ����"H��+������o0�+N���b��f�༻g@:P ��@��t;�kux��3`��Ͱ�$�0�W�ؾ���/��t�h�f�о����<�ؾ� ���7���搾9w�P�ç@�������[� ��� �}оb�ؾT���3H�8L���0��� ��E �(��l0�, �o� ����0�0�H� 8�� ��`���@�30x�N%@�R���\`���@�� �	< �k|�� ����"�� c`���p��L0�z�������.u���:�v �p�v(�+� �'Ѩ;�,��t �����ܾ|����<��tt��	ؾ
��樾<5о�hp������ �����D�(�&!��NȾ��ȾU�Z����kp���侫�Ծ��H�S��`����V��?& ��Ծ�� �y�@�M��O �hc ��[:R� ��� ��8`�S�18�U0��ؽ�x��d��4�S���3� �� �EH���;�V ��H ��� �|I`��p��Ⱦ�!t�H�8�����`�+@`��P�_`�f	`��# �<�`��4��� �|@�����@�w����н�� �gy��Ӣ�����������p���T�(��\ ��{��H�H�z�p��ܾ��ܽK ;����e�ཱི+��2���8Ծ�������e�P��kP��8�J�Ⱦn���o���ճ0�F��6X�[��t����@� 4P�	HȾ^&�g���\�ؽ��0�Q� �䰽�ཥ1`�{��$G �pX��g�:` � ��;��& ������OH�F����p��u�辶�,��u��S(�?�X�_M��X ��}о6XH�"�x��q��V.������2���нӄ �
��x�(��o��E������������Ku �K�x�Wkh�k�辜�<����f�ٳp� {@�@� �  ���$�н <- ����e���K\@����r%о��#��C���~��u,H��'@�\� ��]p��} �=�� r(��ڠ���@:6X �� �9� �=e �8�����0��D��
�о2�0����������P�0vȾf�����@��r0�� ���𽣯0��@�����X0�����, ��<���q����6  �L}��&� �6��@�0��Pp;9 <Bn@���<xh����Y^��P�Ⱦ.����# ��@�(� �ߓ��x^ ��������h� ��耽\�@��R���t�q�@�=Z �p �]>������<bg����p�z`�ȥ ���@����+-����ཌྷMлѮ �1|��v��y�`�� ��  :�< ��=�'� ����⨾W}X��� �&������;��B��J}��� �w�@�]�݈�"���� ��H �B���� �� ��, <�3@<T��7� ����f���fA ��΀��  ��À��0 9� ��H �:6�����
w`���P��逺�$ ��c@�ZY�;�� ��`������UD ��D �g@�f��( ��  9p ��� �L@����"����� ��� ;j �/� ��NP��j��� ;�π;.� �4� �Jh �;࠼��@�j���x���.� :�" �V� �A���bн�2`� ���+�H� =��\ �3����`��� ���P�l� ��^��ex ��� ��`��e ��� �4�`�
:`��w@��P������' �� �:���&`�\8@��0 :� �B�`����t���*p�\$ 9!  �̮@���P�d쀼�� ��� �I.�;M� :  �BG �VW �� ��9���-��r��d �� �� ��� �	ླྀl��:�h�0�#Sར��i@<B�@��� =��� �R� �C�`�J`�ۉ@��〺� <$*@��C �B@�c�@��&����������^h���0��� �!4 �J{ 7�  ��\���̀��� ;�� ;7� ��� �C��f} ;� ��@������4��됼�X@�� �����E7`<"@<�@:�n �*? ��0�;�� ;W( ��h���IP�̆�ޒ�;�U���� �O� �Ƨ���@��� ��< :�p � <8% :�D �Y ��� ��� �h_ �� �[� �54 ���;S� �'����@� �����) ��m�:�� ;�� ��6 �j��;�)��|n �6�оG��G�`��0��<��Z0@�*H �-� �	 ;��<�"��� :� �U��<'*�;����7���)`��^���I�� �{��:�d �>3 ;�5�;�� <�����@��U 8  �&��� �� ;Q� <^�@������@��0��� :�( �B�`��Z�������� �}gv�4� �Q���� �����*�;�Ҁ��B��Z��������R <B���� ���6% ��ʀ��G��� �� �@� �&� ���:�� 6D  :+ <�`�G���^WH���0����+н��0�X�`� 2 �
�`9�� :�h =R@�Y��( ��� ��� :�� �> ��S�*ɠ�*���X�`�H���~0����`���`�T����P��9�/�b���<@�y��:�� 9v� <,� :�
 ���������,��.���5��L~ �M��N�0���ཤ��Rp��M��������U=���p0���`�i����& ��o ��l@�z�CI�;�s ;E0 �p� ;�a �d� �$���8� �ྈ�`�Zq��;ؽ-�`�Y��F� �� ��T��~��bV ;~a �`` �B���]L �@� ��"`��  �����|a�����;�}�9�H �q <� ��� ��
����;׋ �d? ��@�Y����#_X��� �c �Dp��� �>k��D�`��X�H��@�ȾV�����3'��i���q�`�'��N1 ��н�� �'���n����нi*�dv������w���4@� ��  :i� �8� �o� �'�༉	��E` ��u��#轪������������ H�- ��� �kR����ܼ ��� ��� ��� ;rO �ͤ��x������;\� ��� �z�;iY ��n�:)4 ;�� ���z ��r ��1 �F � �Ҡ���� :	 �Qt��c� ������ ���� �8� ���J�R� ����;�j�<�@��j@��&@�.� ��,����D�2=��"� �հ�^���80�1O`��'��@�������#`�,��Z@�� ��Y��F� �� ��� :�� �1H��$� �#H ��ǐ9f@ ��| ��Z�X ��7����E�<3� ;�#���� :UT ;�� <R<@��P �����:����C �ހ�	� <H <ag <�@� �s �����wp������������$ ��`����ƨ �<㠹�� �K~ ���@�ON��S|P�pj`�>π���p��{@��΀�{ ;CX ��D0��n��g ��h���н�н�0��Q@��h�r� ��& �U����p��mཤ� �͞ �6�;b� �8� ��� ��� ��R��� 9�� ��׀�P�@�6��D	 ���`������𽊧P;� ��D��r�@�DR ��p :Ô ��r@;O ���0��� �Z  �F`�- ��Y@�6� :� ��a���� ��3 �j� ��90�8: �PP �5k����@ �v� �YcཛྷL@��`�E�8�\�(���ԾQ�X�$��g�8���>�Ƚ{�@�b �U���`����Q������� �7�H�L�`�HY��8��P�Mh�.��)7`��	���1���7���^ �:J �J� ��� :�" ��f ;� ��X ��  �{U���� �؃���Gн��P�갾q����H�����A�e��� �D4@�� <�� ��@�������P����ݝ��
���F� �{���j@�@耽�	�� ��,. ��� �� �3X�� <Hˀ��z �c� �K?����`��6,�ED��3���[� ��| �U���'�����ʒ��H��kмje���簾:�(�A�X�DE��4�ؾK�P�~�оX���ཪ4 ���p�����2���� ;d �W@ ��N �J��+u��&0��9���� :�� :� �^����'��ٰ�������D� �=�=�x��`��X�'@��K �Q  ;�c ��^ ������Y �����ľ@��(��� ���@<2@�{4 �3\��dc@�X@;u �� ��� ��X�C܀�e�н�Pм9N ���� �� ���K���><���h �W� ��� �T۠����M& �]����� �;� ������� �$q ����I����2〼�*�������d���G@��� <���Ph���8۠ ��O@��\ �o���Rˀ��n ���@��2��6Sཚ���-� ��� ��������<_��H��S���eh� ;Qd ��e��� ��N����������(;����<����  ����@<��`�� <�d�:�J 9r� <8��;�����@��X �0�� �����
`<s�;�ǀ;�K���� ��� �ذ���D��b�9�0 �c2 �A ;� :�b ��� �� �` �J� ��� ����:�@��r��?t �����k �D� ;(� �X�@��8 �#�:�t ;} ��� ��� <��<լ�<?G��$�;�; �� �O������� ��@��� �1' �Sl����`�������@�����&D���� �o����;��:� ����G� �:݌ �_x ��H@��ʀ��4�9�� :� �::`�b	X���Ng�u} ������� ;O; �� �D ����[���P� 9  �n� ;u <8K��� �Eh �8� ��� ;�� 9�� �-� ���н�����:CD ��| �U� :;X :�� �� ��S��C��</_ �m �KZ �� �� ; i �� ��,@����V���~(���B �6���)��D��� ��� ����Sh ��������� �����yP�"x �2�����+�h��5p��,@;�$ <Kb��F� ;{r �+q��J�л� �`�8��~(�+V(<H,@��� ���p������� �0�@������|P� `���ȾI��;Kp�!�@�˿���P�V���l\X���P��0 ���@�� ��� ;C ;���<� :�F �& ��� �e �z �I� �Qx �"� �N �L@ ;����]� �7� ;�C��.� �1D ��T���� �<n ��� ;�< �
$ �]�կ�:> ��� � ��t ��y���� ��4��P�z������;D� �V� ��� � ��4@���`;�� ��� ��, ��� ��@��< <��<�9��:!��T� ��X�L��̴���^�*?���Ⱦ$x������������`��o����8��H�>&��$� ��; �X��� ������� �& 8�� ����9� �� ��值��������A ;�a :� �2H �����} ; d �� �������1�`��^@��u@�.� ��� ��� �0���zƀ��Ā��� ��8 �?@ �� �� ��� ��C��X  <�@�px@�푐��|��#hH�=GH�s� ;� �݃��� ��� <���B}�:&� �f� ��1 �8  �Z,���� :�v �UP �K��9� �M^༲�:� ��p ��Q@�P���	Q �ҹ :)� �5� :�� ;�k ��" �
 <>� :Ch ���@�aP��^ ������ �s> �\� �b  �� ��Հ��� ��# :�V ���@��.@��y �( ��1 ��Z ��h ��� �����ڀ:� ��� ����˼ ��  ��< ��`�w ����� (�������=`���P��𼹏��� �`% ��G ��� ;� ��[�T�����Y� :�f ;OI :�* �� ��70�&Q��b����9��1��� �>~ ;%� ;�� ;�� <D����� ��� ��� �ϟ�������(����p��v��i� ��p�8��;�� �Р��0�G� ��1 �, ����9�p ;����7/ ����� W �� @������8�d- �K���� 7�� ��A ;dg ����� ����q@�y1@��O`��k��ᠾؽ�����p�j ��p��N�:I༒���� �.� �S� �k6��� <{ ��k�;( �$v �s ��4 :4l ��n������F�8�! �ڔ��ݭ 6t  ��"@��@�� �Q| ����>Ƞ��s �����-a��� ��] �S\ �I�@��� �ap ��# ��� ��p0�B����z ��Q ��� ; 9�� �� �c�@��Y������  �F���e(�:G���D0�E���3����=�`<}���Z����(�@㠾��T���`�‽��ՎP��s ��`��T��ˀ��x�q�@�3���PY@����`�(�`:k� �` ;	3 ;j� �r��� �m� ��� �! <�@���:�@ �����ϐ�k����-ལ���n@<,@�,� <��`<���<@@<5/�;�7�����t �_F �4��ڸp�/2����3`�
`�(�:�\ ;�Z��X ;�р:�� �:�@������� ��
 ;Hb �mR ��G`��;�=b�P;y� �3h�<G-?��=�x�Բ��n�оyl�u��8��� ��
 �r`�Or�����L���� �>� �7A�����8 �,� ��o <H�@<nj�;�] 9�� �T :D ��D �y� <�@��� ����Ag��\@ �<;� � ;mq �yR �]Z �N �J���d$��,W ��( ��@�� ��R��c� �r�X�<7�<������!�����Ő�V���1����,����@�� 9i� ��_��Ă �C��;�� �.W��� �Á ��� ��� �� ���;����� �`�;�f �� ��*@���6P�D���a���� �Ғp��h����p�H���V ��@��� �^6����P�?������� �����Ġ:�( 9I� ��� :�@ �Ж@�U� ��!��7���Ŋ@���@9�P ��  �G� ��� ��� ;��;a �� �U� �Đ�[���1 �YĀ���@�C
P������0��x�P(��8�
���| ��a���a@�k�`��� �&@ �%x �O~��) ����9RP ����:� �08���K :\ �@ �; <���<�`:�f ��� �v���n��E�@�QR��V ���˧���"0��v��.����*P�Ό������(������+�������=@�� �z� �>4 �\ ;M� �p  :�< �� 9� ;�.���� �E ;؄ ;�〸�� :�> ;�~ �*쀼R���  ��p�� ���@;������ ��Oн��@�޺��ü�=�p �|@���Y� �Ǖ ����6���$2 ;F� 9�� ��y ������������������/𼯤��/���)%@�\�`�����]/ ��&@������KP�h7�ҍ��ň�t| �Љ�;ۿ <G@�� �W;�U-�y�����0�m2���jp��D ��q ��2`��� �is�9�  �4%Ƚ�6�<_� �=���,�
� �����o= ����������������=6b���@��� �π ����<� �R  =8&����@��������� =��8�� 9F� :_� �:� �
ep����;�| ;刀;�=�����a@�������������������������`;�a �Z& �d���� �oo �k�+ ��e@<9E��]F���� :�� ;�� ��оHv �0��=U�0��
���������)���D=���>�|=��ȿ�h�������������������������nVP���4�����������������������������������������������������������������������������������������Q=�P<u^�=]@����������� ��� ���������������������������������<̐����@���=����*<p) ����D��V�p��{`���l���𾳃������Gɾ������p�ћ�=tBp=�P(=��>"oh�����"�X�+�X�A�������!�=բ�>18��1�>>vT=��<p"�����������s�ׁ =pY�������PV���������������������������������������p������w0��� �I�༇���)K���Ȅ��u\����������Ub����
V���������������������������ཱི�����