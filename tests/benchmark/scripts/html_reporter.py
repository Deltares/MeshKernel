import os

from airium import Airium


def generate(work_dir, plotter_report):
    page = Airium()
    page("<!DOCTYPE html>")
    with page.html(lang="en"):
        with page.head():
            page.meta(charset="utf-8")
            page.title(_t="MeshKernel Benchmark")

        with page.body():
            # with page.h3(id="id23409231", klass="main_header"):
            #    page("Hello World.")

            for graphic_path in plotter_report:
                with page.div():
                    page.img(
                        src=graphic_path,
                        alt="alt text",
                    )

    html = str(page)  # casting to string extracts the value

    print(html)

    with open(os.path.join(work_dir, "report.html"), "wb") as f:
        f.write(bytes(page))
