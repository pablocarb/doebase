from unittest import TestCase


class Test_import(TestCase):
               
    def test_import(self):
        b = True
        try:
            from doebase import(
                makeDoeOptDes,
                callDoE,
                evaldes,
                doeRequest,
                getDoe,
                mainDoe
            )
        except Exception as e:
            b = False
            print(e)
        self.assertTrue(b)