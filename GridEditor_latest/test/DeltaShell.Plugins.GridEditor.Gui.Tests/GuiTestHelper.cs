using System;
using System.Drawing;
using System.Drawing.Imaging;
using System.Runtime.InteropServices;

namespace DeltaShell.Plugins.GridEditor.Gui.Tests
{
    public static class GuiTestHelper
    {
#if Generate
        private static bool generate = true;
#else
        private static bool generate = false;
#endif

        [DllImport("msvcrt.dll")]
        private static extern int memcmp(IntPtr b1, IntPtr b2, long count);

        public static bool CompareBitmapWithReferenceFile(Bitmap b1, string referenceFilePath)
        {
            if (generate)
            {
                b1.Save(referenceFilePath);
                return true;
            }

            using (var b2 = (Bitmap) Image.FromFile(referenceFilePath))
            {
                return CompareBitmaps(b1, b2);
            }
        }

        private static bool CompareBitmaps(Bitmap b1, Bitmap b2)
        {
            if ((b1 == null) != (b2 == null)) return false;
            if (b1.Size != b2.Size) return false;

            var bd1 = b1.LockBits(new Rectangle(new Point(0, 0), b1.Size), ImageLockMode.ReadOnly, PixelFormat.Format32bppArgb);
            var bd2 = b2.LockBits(new Rectangle(new Point(0, 0), b2.Size), ImageLockMode.ReadOnly, PixelFormat.Format32bppArgb);

            try
            {
                IntPtr bd1scan0 = bd1.Scan0;
                IntPtr bd2scan0 = bd2.Scan0;

                int stride = bd1.Stride;
                int len = stride * b1.Height;

                return memcmp(bd1scan0, bd2scan0, len) == 0;
            }
            finally
            {
                b1.UnlockBits(bd1);
                b2.UnlockBits(bd2);
            }
        }
    }
}
