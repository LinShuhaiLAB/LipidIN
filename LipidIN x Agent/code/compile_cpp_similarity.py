import os
import sys
import shutil
import tempfile
import subprocess

MODULE_NAME = "eq3_similarity"


def compile_similarity_cpp(cpp_source_path="eq3_similarity.cpp"):
    """把 eq3_similarity.cpp 编译成可导入的扩展模块，编译产物复制回当前目录。"""
    print("=" * 60)
    print("C++ 相似度模块编译")
    print("=" * 60)

    cpp_source_path = os.path.abspath(cpp_source_path)

    if not os.path.exists(cpp_source_path):
        print(f"✗ 找不到 C++ 源文件 {cpp_source_path}")
        return False

    print(f"✓ 找到 C++ 源文件 {cpp_source_path}")

    original_cwd = os.getcwd()
    cpp_filename = os.path.basename(cpp_source_path)

    with tempfile.TemporaryDirectory(prefix="cpp_build_similarity_") as temp_dir:
        print(f"✓ 创建临时编译目录 {temp_dir}")

        shutil.copy2(cpp_source_path, os.path.join(temp_dir, cpp_filename))

        setup_content = f'''
from pybind11.setup_helpers import Pybind11Extension, build_ext
import pybind11
from setuptools import setup
import platform

if platform.system() == "Windows":
    extra_compile_args = ["/EHsc", "/bigobj", "/O2", "/std:c++17", "/utf-8"]
else:
    extra_compile_args = ["-O3", "-std=c++17"]

ext_modules = [
    Pybind11Extension(
        "{MODULE_NAME}",
        ["{cpp_filename}"],
        include_dirs=[pybind11.get_include()],
        extra_compile_args=extra_compile_args,
        cxx_std=17,
    ),
]

setup(
    name="{MODULE_NAME}",
    ext_modules=ext_modules,
    cmdclass={{"build_ext": build_ext}},
    zip_safe=False,
    python_requires=">=3.6",
)
'''

        with open(os.path.join(temp_dir, "setup.py"), "w", encoding="utf-8") as f:
            f.write(setup_content)

        os.chdir(temp_dir)

        try:
            subprocess.run(
                [sys.executable, "-m", "pip", "install", "pybind11"],
                check=True,
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
            )
            print("✓ pybind11 就绪")
            print("开始编译 C++ 相似度模块")

            subprocess.run(
                [sys.executable, "setup.py", "build_ext", "--inplace"],
                check=True,
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
            )

            print("✓ C++ 相似度模块编译成功")

            compiled_files = [
                f for f in os.listdir(".")
                if f.startswith(MODULE_NAME) and f.endswith((".pyd", ".so", ".dll"))
            ]

            if not compiled_files:
                print("✗ 未找到编译生成的文件")
                return False

            for compiled_file in compiled_files:
                shutil.copy2(compiled_file, os.path.join(original_cwd, compiled_file))
                print(f"✓ 复制编译文件 {compiled_file} -> {original_cwd}")

            print("=" * 60)
            print("🎉 C++ 相似度模块编译完成")
            return True

        except subprocess.CalledProcessError as e:
            print("✗ 编译失败")
            print(e.stdout)
            print(e.stderr)
            return False

        finally:
            os.chdir(original_cwd)


if __name__ == "__main__":
    if compile_similarity_cpp("eq3_similarity.cpp"):
        import eq3_similarity
        print("✓ 模块导入成功")
