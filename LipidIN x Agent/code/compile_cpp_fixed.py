import os
import sys
import subprocess
import tempfile
import shutil
def compile_with_fixed_path():
    """Compile EQ3.cpp as eq3_search."""
    print("=" * 60)
    print("C++模块编译检查 | C++ module compilation check")
    print("=" * 60)

    import os
    import sys
    import re
    import shutil
    import tempfile
    import subprocess

    module_name = "eq3_search"
    cpp_filename = "EQ3.cpp"

    current_dir = os.getcwd()

    existing_pyd_files = [
        f for f in os.listdir(current_dir)
        if f.startswith(module_name) and f.endswith((".pyd", ".so"))
    ]

    if existing_pyd_files:
        print(f"✓ 检测到已有编译文件 | Existing compiled file detected  {', '.join(existing_pyd_files)}")
        try:
            import eq3_search
            print("✓ 已有 eq3_search 模块可正常导入，跳过编译 | Existing eq3_search module can be imported successfully, compilation skipped")
            return True
        except ImportError as e:
            print(f"✗ 已有编译文件存在，但模块导入失败 | Existing compiled file found, but module import failed  {e}")
            print("继续重新编译 C++ 模块 | Continue recompiling the C++ module")

    else:
        print("未检测到已有 eq3_search 编译文件，开始编译 | No existing eq3_search compiled file detected, compilation started")

    cpp_source_path = r"D:\bio_inf\02 LipidIN 3.0\pycode\EQ3.cpp"

    if not os.path.exists(cpp_source_path):
        print("✗ 错误，找不到 C++ 源文件 | Error, C++ source file not found")
        print(f"  路径 | Path  {cpp_source_path}")
        return False

    print(f"✓ 找到 C++ 源文件 | C++ source file found  {cpp_source_path}")

    with tempfile.TemporaryDirectory(prefix="cpp_build_eq3_") as temp_dir:
        print(f"✓ 创建临时编译目录 | Temporary build directory created  {temp_dir}")

        try:
            temp_cpp_path = os.path.join(temp_dir, cpp_filename)
            temp_setup_path = os.path.join(temp_dir, "setup.py")

            with open(cpp_source_path, "r", encoding="utf-8") as f:
                cpp_content = f.read()

            cpp_content, replace_count = re.subn(
                r"PYBIND11_MODULE\s*\(\s*[A-Za-z_][A-Za-z0-9_]*\s*,\s*m\s*\)",
                f"PYBIND11_MODULE({module_name}, m)",
                cpp_content,
                count=1
            )

            if replace_count != 1:
                print("✗ 未找到唯一的 PYBIND11_MODULE 定义 | PYBIND11_MODULE definition was not uniquely found")
                return False

            with open(temp_cpp_path, "w", encoding="utf-8") as f:
                f.write(cpp_content)

            setup_content = f'''
from pybind11.setup_helpers import Pybind11Extension, build_ext
import pybind11
from setuptools import setup
import platform

if platform.system() == "Windows":
    extra_compile_args = ["/EHsc", "/bigobj", "/O2", "/std:c++17", "/utf-8"]
    extra_link_args = []
else:
    extra_compile_args = ["-O3", "-std=c++17"]
    extra_link_args = []

ext_modules = [
    Pybind11Extension(
        "{module_name}",
        ["{cpp_filename}"],
        include_dirs=[pybind11.get_include()],
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
        cxx_std=17,
    ),
]

setup(
    name="{module_name}",
    ext_modules=ext_modules,
    cmdclass={{"build_ext": build_ext}},
    zip_safe=False,
    python_requires=">=3.6",
)
'''

            with open(temp_setup_path, "w", encoding="utf-8") as f:
                f.write(setup_content)

            print("✓ 文件准备完成 | Files prepared")

            original_cwd = os.getcwd()
            os.chdir(temp_dir)

            try:
                print("安装 pybind11 | Installing pybind11")
                subprocess.run(
                    [sys.executable, "-m", "pip", "install", "pybind11"],
                    check=True,
                    capture_output=True,
                    text=True
                )

                print("✓ pybind11 安装成功 | pybind11 installed successfully")
                print("开始编译 C++ 模块 | Start compiling the C++ module")

                result = subprocess.run(
                    [sys.executable, "setup.py", "build_ext", "--inplace"],
                    capture_output=True,
                    text=True,
                    check=True
                )

                print("✓ C++ 模块编译成功 | C++ module compiled successfully")

                compiled_files = [
                    f for f in os.listdir(".")
                    if f.startswith(module_name) and f.endswith((".pyd", ".so", ".dll"))
                ]

                if compiled_files:
                    for compiled_file in compiled_files:
                        dest_path = os.path.join(original_cwd, compiled_file)
                        shutil.copy2(compiled_file, dest_path)
                        print(f"✓ 复制编译文件 | Compiled file copied  {compiled_file} -> {dest_path}")

                    print("=" * 60)
                    print("🎉 C++ 模块编译完成 | C++ module compilation completed")
                    print(f"生成的文件 | Generated files  {', '.join(compiled_files)}")
                    return True

                print("✗ 未找到编译生成的文件 | No compiled output file was found")
                print("编译输出 | Compilation stdout")
                print(result.stdout)

                if result.stderr:
                    print("编译错误 | Compilation stderr")
                    print(result.stderr)

                return False

            except subprocess.CalledProcessError as e:
                print("✗ 编译失败 | Compilation failed")
                print("标准输出 | Standard output")
                print(e.stdout)
                print("错误输出 | Error output")
                print(e.stderr)
                return False

            finally:
                os.chdir(original_cwd)

        except Exception as e:
            print(f"✗ 编译过程出错 | Error occurred during compilation  {e}")
            return False
if __name__ == "__main__":
    success = compile_with_fixed_path()
    if success:
        print("\n测试导入模块...")
        try:
            import eq3_search
            print("✓ 模块导入成功!")
        except ImportError as e:
            print(f"✗ 模块导入失败: {e}")
