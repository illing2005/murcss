CDF       
      lon    H   lat    $         CDI       GClimate Data Interface version 1.6.1 (http://code.zmaw.de/projects/cdi)    Conventions       CF-1.4     history      �Wed Jul 09 12:03:24 2014: cdo sub /tmp/murcss/cache/constantField.nc /tmp/murcss/cache/2_5_tas_miklip_uninitialized_mpi-esm-lr_uninitialized_1960-2000rmss_sqrt.nc /tmp/murcss/output/20140709-120305_tas_miklip_initialized_mpi-m_mpi-esm-lr_initialized_miklip_uninitialized_mpi-m_mpi-esm-lr_uninitialized_HadCRUT_1960-2000/2-5/miklip_uninitialized_mpi-esm-lr_uninitialized_input2/msss/2_5_tas_miklip_uninitialized_mpi-esm-lr_uninitialized_1960-2000_rmss.nc
Wed Jul 09 12:03:24 2014: cdo sqrt -sub /tmp/murcss/cache/constantField.nc /tmp/murcss/output/20140709-120305_tas_miklip_initialized_mpi-m_mpi-esm-lr_initialized_miklip_uninitialized_mpi-m_mpi-esm-lr_uninitialized_HadCRUT_1960-2000/2-5/miklip_uninitialized_mpi-esm-lr_uninitialized_input2/msss/2_5_tas_miklip_uninitialized_mpi-esm-lr_uninitialized_1960-2000_msss.nc /tmp/murcss/cache/2_5_tas_miklip_uninitialized_mpi-esm-lr_uninitialized_1960-2000rmss_sqrt.nc
Wed Jul 09 12:03:23 2014: cdo sub -sqr /tmp/murcss/output/20140709-120305_tas_miklip_initialized_mpi-m_mpi-esm-lr_initialized_miklip_uninitialized_mpi-m_mpi-esm-lr_uninitialized_HadCRUT_1960-2000/2-5/miklip_uninitialized_mpi-esm-lr_uninitialized_input2/msss/2_5_tas_miklip_uninitialized_mpi-esm-lr_uninitialized_1960-2000_correlation.nc -sqr /tmp/murcss/output/20140709-120305_tas_miklip_initialized_mpi-m_mpi-esm-lr_initialized_miklip_uninitialized_mpi-m_mpi-esm-lr_uninitialized_HadCRUT_1960-2000/2-5/miklip_uninitialized_mpi-esm-lr_uninitialized_input2/msss/2_5_tas_miklip_uninitialized_mpi-esm-lr_uninitialized_1960-2000_conditional_bias.nc /tmp/murcss/output/20140709-120305_tas_miklip_initialized_mpi-m_mpi-esm-lr_initialized_miklip_uninitialized_mpi-m_mpi-esm-lr_uninitialized_HadCRUT_1960-2000/2-5/miklip_uninitialized_mpi-esm-lr_uninitialized_input2/msss/2_5_tas_miklip_uninitialized_mpi-esm-lr_uninitialized_1960-2000_msss.nc
Wed Jul 09 12:03:23 2014: cdo timcor /tmp/murcss/cache/tas_Amon_mpi-esm-lr_uninitialized_r1i1p1_196001-200512.nc_1961-1965uninitializedBQZ4MK.remapped_ENSMEAN_ANOMALIE_selYear_2-5279960Y2YIZTimmean2_5_tas_miklip_uninitialized_mpi-esm-lr_uninitialized_1960-2000AEVHGA_merged.nc /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1961-1965Q3A2P4.remapped_ANOMALIE_selYear_2-5277605NLS7LTimmean2_5_tas_miklip_uninitialized_mpi-esm-lr_uninitialized_1960-20002ZA9R8_merged.nc /tmp/murcss/output/20140709-120305_tas_miklip_initialized_mpi-m_mpi-esm-lr_initialized_miklip_uninitialized_mpi-m_mpi-esm-lr_uninitialized_HadCRUT_1960-2000/2-5/miklip_uninitialized_mpi-esm-lr_uninitialized_input2/msss/2_5_tas_miklip_uninitialized_mpi-esm-lr_uninitialized_1960-2000_correlation.nc
Wed Jul 09 12:03:23 2014: cdo mergetime /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1961-1965Q3A2P4.remapped_ANOMALIE_selYear_2-5277605NLS7LTimmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1966-1970CH76ZC.remapped_ANOMALIE_selYear_2-527761JUPMVETimmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1971-1975C6WG6O.remapped_ANOMALIE_selYear_2-5277623HN4WDTimmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1976-1980ZGP54Q.remapped_ANOMALIE_selYear_2-5277637Z7LTWTimmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1981-1985T74HT8.remapped_ANOMALIE_selYear_2-527763AIZEZ7Timmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1986-19902WG1Q7.remapped_ANOMALIE_selYear_2-5277601KLFUFTimmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1991-1995IQC7AZ.remapped_ANOMALIE_selYear_2-5277624EBRR2Timmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1996-2000G2YF3T.remapped_ANOMALIE_selYear_2-52776172OAT1Timmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_2001-20053Q7AVD.remapped_ANOMALIE_selYear_2-527763YQWBFXTimmean
Wed Jul 09 12:03:22 2014: cdo timmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_2001-20053Q7AVD.remapped_ANOMALIE_selYear_2-527763YQWBFX /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_2001-20053Q7AVD.remapped_ANOMALIE_selYear_2-527763YQWBFXTimmean
Wed Jul 09 12:03:22 2014: cdo seltimestep,2,3,4,5 /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_2001-20053Q7AVD.remapped_ANOMALIE27763yearmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_2001-20053Q7AVD.remapped_ANOMALIE_selYear_2-527763YQWBFX
Wed Jul 09 12:03:22 2014: cdo yearmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_2001-20053Q7AVD.remapped_ANOMALIE /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_2001-20053Q7AVD.remapped_ANOMALIE27763yearmean
Wed Jul 09 12:03:19 2014: cdo sub /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_2001-20053Q7AVD.remapped /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1986-19902WG1Q7.remapped_2000_obsYT1ZHG /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_2001-20053Q7AVD.remapped_ANOMALIE
Wed Jul 09 12:03:19 2014: cdo ensmean /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1986-19902WG1Q7.remapped /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1991-1995IQC7AZ.remapped /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1961-1965Q3A2P4.remapped /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1996-2000G2YF3T.remapped /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1966-1970CH76ZC.remapped /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1971-1975C6WG6O.remapped /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1976-1980ZGP54Q.remapped /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1981-1985T74HT8.remapped /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1986-19902WG1Q7.remapped_2000_obsYT1ZHG
Wed Jul 09 12:03:12 2014: cdo remapcon,r72x36 /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1981-1985 /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1981-1985T74HT8.remapped
Wed Jul 09 12:03:06 2014: cdo selyear,1981,1982,1983,1984,1985 /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc_1981-1985
Wed Jul 09 12:03:06 2014: cdo chvar,temperature_anomaly,tas /home/illing/workspace/murcss/sample_data/HadCRUT.4.2.0.0.median.nc /tmp/murcss/cache/HadCRUT.4.2.0.0.median.nc
Updated at 26/06/2014 15:02:28   CDO       HClimate Data Operators version 1.6.1 (http://code.zmaw.de/projects/cdo)          lon                 standard_name         	longitude      	long_name         	longitude      units         degrees_east   axis      X        @  h   lat                standard_name         latitude   	long_name         latitude   units         degrees_north      axis      Y           �   const                      
_FillValue        ����   missing_value         ����     (�  �        @      @$      @.      @4      @9      @>      @A�     @D      @F�     @I      @K�     @N      @P@     @Q�     @R�     @T      @U@     @V�     @W�     @Y      @Z@     @[�     @\�     @^      @_@     @`@     @`�     @a�     @b      @b�     @c`     @d      @d�     @e@     @e�     @f�     @g      @g�     @h`     @i      @i�     @j@     @j�     @k�     @l      @l�     @m`     @n      @n�     @o@     @o�     @p@     @p�     @p�     @q0     @q�     @q�     @r      @rp     @r�     @s     @s`     @s�     @t      @tP     @t�     @t�     @u@     @u�     @u�     @v0     �U�     �T�     �S`     �R      �P�     �O@     �L�     �J@     �G�     �E@     �B�     �@@     �;�     �6�     �1�     �)      �      �      @      @      @)      @1�     @6�     @;�     @@@     @B�     @E@     @G�     @J@     @L�     @O@     @P�     @R      @S`     @T�     @U�     ��D8��g4�����������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������������KE��� �������������������������������������hǘ��0��Ŀ=�L�����؈��O��\��"���a�+TV�=* ����������������=�k>1�����������������������������������������>�F�>�庾dx=�����B����DL��D�?�@������������=������������L����������������������������������������������������������������������������������������������������������������=��>"�\�0����P�E��|/���@�%U�V	v��yH��,�s,�V��𐽺�p��䤿D����Z���̈��䂿/K$��o���W@��� ���P������������>�#�������y��Y&`�(Z���7\���ʿ������8��ԑ �+(���t���������Lx��)�`����H^H���侼P�H�0�LX��`�(�;���x��P � X>H=Qπ=�����̾,�5`>+c�= ��=�kX>������&�j>�p~>S!��]#8�:Z��~J���F���Ŀ�FD����������iH�D,οQŨ�=־�/�������\��p�q~�`���4�����'�޿~�¿oʲ���񙴾�I����4�*�P�8>I�=R_ ���T���佭���,�����?��J�P��a��7Ø��'��)�X�'�
��*����d�碨��uԾOr(�Gp� ؽȀ�6,(���Ⱦ�͜�b���н���d�8����)�����%�8���������� �  о���}���}G辢�3�����4�^ਾ��X�2rؾ$��d� �t�]����݀�P >(� =:p�����X��oH�z� �ɓ <S�������߼������4��]P��:����H�� �o�8������D�T�о��d�F�D����B���
ા��Ծ��(�Q���ͤ<�H���м�J ��� ��� ��0�:+ ��ڀ������C`��&���W0��b������I�H�~X =����/Q`�� <&� ������@�16P�*UP�o� =}� ;�X��pW��I�<��=�<���<���f��N9x�l�о$����оP�X�У�:�@��܀�[��Il �_� ��p���@=���<4�@��߀��� �?
�=�r=}�=���=���=c�=OP��W ��m���L�Q����踿�C���6��-*���@�/�>��a4�O�v��s�;�} <:������M����8 �� =��X� ���9@��P�n����$<�`>�0<�$ ;�% <����$ <�� =*�p=��P=�	��J�@��0=Y�@>�T>P�\=_p�=�0`>"�=3/ =�@>BH����bQ0�+V�V�����>�(>ᠾ�`��	���9��������=倽�|���D���H�V�`�Jf ��#�;s� �|� �#� ��� =���>IFx>�!�>���>T��=�F ��0�@%п��Z���F�t <'���d =���=�3=26 �|䰾&�P��ǰ�`�i����@<� =��h>p>PwT>�-�>Ԯ�>lt,>@��=���>9t�H��=�}h=�x=�� =�0P��ཨ&0�BȀ=��H�P �6툾�Ԩ�%�>�JX���ھ�J���]`�Cy@��X��+n�bv���BX�;�X���=�X=�6 >NH=Q�`9i  ��  ��� �W �m��"␽�J�=�`=���=���<�A���A ��:�=GOP=���<J�>���>��>.� >[=��<�,����<�Ӏ�X����� ��# ��hP�Ꞿ��p���|�8�X�CoB> \|>�H�>��>u),>"V�>�`>�5d>(�>��h??}>���>, �>e�p>Rg(�J�@���p��𾖦4���@� �y��D� ���=���=m��>�7P>��= k�����=��N
�a�>/�>̰=�t(=�?�=��(=E�`=�����{��
@ �����ź�
Q��d@=��>)`=�Z>���� ��?��X�X����7z���2>��2>Mw�>Q�>
�T�4 >^c�>�2�>�������ڰ�7y���} ��t�>�>�*�>���>�m����=���>E|=��>���>�^>��D>�[�? D�>��>��>��P>��6>O\ܻ� �V� ����3���G(���̾�#p����G�=�N=�"�>��(=�/=��8>?0H=,Ǡ�F� < o ���@�՞���mP���p�= ��n@���@�ѧ ��6��㌴�����nnH�� ����;N� =�(�60�<H1@�D��Q�������þh������N�<ڌ �c����� =�n�=��0��e�">j@>�2R��Ҡ>���>�!�>�\�>���?"`]>N�\=��>��r>��>�'^>I��>�sd>��F>�q*=�A0>y�>��=\�0<�8 �U� =�V8>�̾:���:P�N02�*����T�;���'0�Ƭ�=��0=�r�>�n>9Rо(�Ⱦw�8�@���p�CP����@<у`=�b�<��`�Ri��ۨ��/0�����<��Ѱ��|о1J@� t>��>�F�=�İ��ۀ�o��
���U`�)��b���L��rȽ�٠��@��L =�b�>:� >7H�>Z�@>��n=��>��>N>�=^�0>uS��� <��`>�	>��>��P>�C�>b+`=�@���`=�J�=�x<�i�=t�"����`<�\@���@> ��=���;��#��=�U�>Q�>���>��ؽ���>z$>;�=���� @�����/�0���\3�><ؾ+�軗3 >Ӡ>'֜=���=a`����?p=9 �aR@�H��0��T��ؾ��4�)P�� Ⱦ]U���_v���Կ~��XRh�.���(��=�R >E�(>a�l>l��>^D�>	P4�z =���=ư ���=G>`>v�=��>���>�=B�>��d>i��>*�=�>��Ā��p�.����d��.�=����3 >N�=o`>o����@����>��<ꁀ=���>J�>d�>4?=�H�Sۀ��a��sH�|;��������s��ҬP=Kb =�<�=&�=�̈=�`>��,>�>�b.>5h$>h�|�7ʠ�b�辰t��l��1о�%�+ؾFm��>���!̽����u�=�`�>���>��>1��>��=�0<o��<�@=<�`>N��+��_����$>��F>�(p>�	J=[J0>�x�i��G�=��>x%�>O�>W��>�[>��>��\>I�$=��>0t�=�������ؼ�@�U
ľ�!��ߞP�?�@��x�o�����>:��=�����K��6�k`�c��<
 ����<�� >;�>�4�>��P>���>jgL>0-L>�X>��=c����྅E$�t+�w@�琐�x�Ƚx�`�E�������̾��=�qH?�?
�?hx=�޸�<��g��=�^ؽ�Z��2@>���� �����y���=)
@>�AF>��>��4>�:,�xm >��>i|=�'P>���>��>ұ:>�>�I�>�"�>��\>׃�>�|l>%@>���>���K����,�>rk�>G���It�X> �L>H�t�<� �@�<�u�;/� =I�p=�X�� =D��MǠ���=�`=�8>��>K�>"=&.�=`ݠ=r�<��<¢��ec`�� �	��Y?����r)�DHڼq >���>�kR>�j*?}�>V^�;�� �J�ྦ����qH���@�u4�>�4>��&>=�=���ʂP���>:�`?�|>�.V���>��L>ިN>���?c�? ��>�f>���?ծ>��"?&�?62>�x?"�?(Q�?$h�?s�?��?H>�b�=�<=�,��n�=��=�f(=�S�>j�>�e0�RP�{8���0�'� �E �/��E ���4 =�� >"�=���= =�-�>	3 =�P=Nཷ?p��o >b���@���P��t�������>5�>_�p>P,>�J�>N�����d��.�r��=�iؽ� >��>޳d>���>��b���4����>��^?��>��"�"����>�=������>�R^?�?�'?�*?
 >ٞ\>J��>�B�?�n?"�s>�߆=v <�(�>�i��/� ���=q���W���W���� �	G0����*K �-���J��C��;q��Nj���D��F��0�ZG����d��(����T�;���FG�=��=�� >$�>C��>Oc�>7kh>Q���Q�h�)��=4�=�ɀ>�	�>���>�k>��4>�Ĕ>Ѵ�=�<�T >~P=��>���>U�t�MͰ�[�x���d����=ҷX? ʸ>�V�����YT�kк���@���>��>�L?�`?��?#0�?�?�"?%��?g<?
Q>�l*�s��Nd> %=�+м�n@��9�?J����L�`9X�Y�8�
�h����������uIx���(�u�о����qԾ�����[�.��2�+�ʿO����e���>h� >�p�>���>�U��@*��1�=��<�H�=?��� >y��>.� >j,>�e"?e>���>�L�>�4>��N>���>�_̼׀�	�%��hp�|U@>�DF>��>�o�=���7Z��I��n�|q���T>��8>��2>�m�>�K�>���>�/�>�s�?�5>ݍ(?&�>���>�X�>��>�>i��>���>�h��Un��MvX�ޤ\�A7 ��h==&@�B@�<�@��F���7l�4�оv˸�-b��z�I����|�k:��h�E������`������_�>|2T>��>�Q�>r��B�H��V��sP=�~��� =۟x>�7>�D�>�/�>��<>�?5�o?�f>�=�>��=�)��= =�c0>���>��`>��8>��8��6������v��[�p�㛬�G8����>��"?r	>���>D|(>�A>�� >�[F>��x>�|�>�db=�y�<a逾�᤽t��>��>3Y$�[,���ԿJ��Y$��{t��(���|���@�Ct�t���l8��n$�3ˠ��(��Yоc���V�H��5,�\޿��h��=��v�����<	���{j�>^k >�[ >EX>�$>1Q�=���=UzP=�e>	�><�X>��D>y�@>p�?
; >�W>n��>�>�f>�T��C<>��l����?��>�y,��@��`>6�d>H!�=��x>�p�� �ڒ�=��P>n��>�����l��8�� �F>l��=���=R��h� �݄ ��C�y�>`v@���ܿ7���^p�~㬿k?��;Iп���@��=t��.H�ʾ��������{h���̾������X��A���1��`L�}3�H ���(������\��9�H<���=��0��m�=�)@>��;]9 �K��=o4�>�>��>��>[�P=��0>�֦>���>��j>��>��>5�|����>��>�o(>�o�>)q�����>�jH>����似� >�>��t��B4��O���<��p>�#�>�辈PH�+X>��>�ˊ��f@�"T���bp�6����>��$�p��P��p�8� �v�¿�j��ɔ���0��U濈��j7��t��W��L��=�z�0a���`�?҈�K����@>J��>87�?��>�
&>�W��#�@��4����������J>G��>���z)��� =���=���>��>�Z>Z�x=�~X>�X����>�8�>т�>���=��@>�9�(�x������|t=���=����{1�>>:,>�p�=�H�=Zk@=��lN�>pr�>=t<?�?#��>���=�F�>��X>ɼ>dF(=��65 <@4�>G�=� �pm������vԿI"���(ں������8��Ҭ��� �	���m�P�k�B������&�^sľ����'��݌��� �����m�>�!�>=  ���=v��>{�X>�=�J�>a��>rۀ=&�н ���L6x���0>T��>��>��=հ>���>��\>�7l>�\�>�К>56=x%�=�=EE =��>��f>���>��x>�#�>��x>Y���) ���8�p?�>���>�P>h��>��0>��F>`�\>@ڼ>i >��><�X>�>~��>���>�;�>���>��v=�U�=�� ��� ���8��ྸ����8�h������7�z�����v��K��!(?r�>��z��������վ��v�P��HT��u,���нߥ�>��8>�_f>�y>�,�>�6�=���H������4=߉>0�>(�D=���>��>>�L>�:>�����T@�x%h�����2l���D=Z#p>F�>���?��?֬>�ް>�>�b�>Ũ�?׿?"Kj>ۊ�>�2�>�;�>PU�;�!�=: �Y ���0���̿i-����T�5~�eJ�ê,���@=ř8=�Q�<�����оr����<���!(��O<��� � ��H��.` >�?"O'?8v>syо�V$�&��R6�1☽�g@>;�<��@=mq@>G��>�߄>�BZ�3����K�: 8���<� @=��=�L�=�P=Hn�>�G�>��.>�0�>���>vut��uȿS���̘��4�=������>��\>�MV>�[�?��>�u">�n�>�B�>��6=�o�>�P^>�נ>��>��|>�S�>��P=yA�t� ��G��m��K����/f��1d�����^!���ꈾ��ؾn�_依k����H������6��H�7/`�c �?�>$>��>�8>�#�>�l���]��x�������,�4(=*��<�	`<  �JK��=񠾍&���Bྐ~��D���=r�0>y��>�k*>���>�1>�x>���>�G$>̯�>��>���>88�P����i�>�>G��>�W�>�%p>N2���<�@>D1,=�H>$� >r��>�4=��P=G��=�^�>��>d6�>C��=Rh�<�Y��0���ݲ���X���t��� ����5h�!2��/?��8�,V�.�ξ�x(�Px>"��>�?H>�B>\[�>W�\>d��>iȴ>,i�=����p�IȾ\������������;�������=�[�>;Ot=�M�jp�v��'ߠ=
:�>� �>��>��>Ώ�>�A�>�5>ϖ�>�b>��`>�s�>�w6>�1>z��>)�>Z0=���<^`���H��7� =�� =۝x>u�H>�k\>��j>�L�>5Ô>w��>�>|q|>��>�,H>�~���
��=�$�B���P�r�̿�I�� �0�PֿY���c	�J4��������鄾y�(�[��=tB�=��P>K|>Q�>9��=��`��b����F���C��Z#����	q�>2�P��tܾ���������;�l �������ݝh����waP��|>�#6>��F>�l�>��>�N>�E�>x�D>k�4>f0�>ch>��>g��>j.x>/�=� ȼMƀ�;W`="�`=����(�<�| >�0><L >Z�|>��>��">��>��>�
d>��$>����w�f���¾�gL����ۈ��jZ��r��+˄�����$��"@��Ⱦ����0�$�<���>��>���>Mc�="���
�B����#t��%p��%�������5��o|��*](��f�)XP=�H���90���̾p=�H��Y� �tRx�]Pp���>p��>멢>��>�v&>�j�>�|>��>>�8>SU�>O�=�~`<s(��TF�ܲ�)���9R��M
��EؾvH�6E�=;��>	8�>M�$>�"�>:��<j��=Rm�>]7>ܞ�>���>�/(>�A�>���>{l�����Y���N�a����r�X;��o���7��l��׎���/p�+��=�4H>*�<>��Z>��=�k�=��@�N��@�����x���.�����|�^
X�!���� ����`=�p��Ⱦc�x���Ծ����4���@轫�@=�G=�d�>1��=�)�>��t��>��?��>֋Z>|U<=�H��o �2������R$���<�X/� 8�] ��A ����dy�;�d =��X� ��%P�'A ��p����=��`��� >��|?�r?�0?B�>���>Z� =FjP<���.(������D��Qܾtx�3ߠ���P� ��>}�>�� >���>�o�>w�>1��=;��9+��A�h�����@=���>s�=5�`>]�8>lQ=�Hн� ��@������Ԁ����>� ��k���H�==�+���о�9�����;� >[�:�< ���$�
��w����]��ƶ��2޴�0-��$��"`<�� <�� �p�������������{��u��%�.�.J�-J0��Qt�5�=���>-R�>x�$>B�>���>�
�>���>w@>!��>Gބ=�" ���@�fl >�.>>���>��>�1T>��>�JZ>��>�-�>�gL>K >��ؾi�P�/��>��(>�Md�( >�e4?��>��>uj�=k@<�����>�����������������>��6>a��7C=>�p>�d��+�Ҿ��==d���(��P�8�`�|�����*��
��%��`����`��j�7��˅X�X����g{�>T��� �p�>i]�>T]��z�,P�t|�>�7&>������>4ZĽ��������H��AӰ��k2��C$��!,�����􄐿��[n0��K0��g,��q������������>�>6>������>%>Q���. �����V0�'b��vt�>Rcx>;G�K��>B�t������������������������>�cl=������H�R��[*@� |N��o޿����,�χ��?z�� �5�4�
�*.��.��=�p="pп�f�?L`�>�Y�>�n^��� �o��+J�>�]R>�2f��������������������������������������ؼ>����(����������Jh�������������v-���b���N��BL�с0����������W<lw�������Qp��j���>���>�	�>�xB<!�@����������谾�r�����������������������������>��$>�H>��&���@���������������������$K�����J��0��Ah���L�����:t� �`�������?�����[��������������������������������������<U:��ؤ�T=�������������t�6�?bJ�����������������������`��ٰ���������������������������������������������������������������������������H����v��������������������������������