import argparse
import importlib.machinery
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
CODE_DIR = Path(__file__).resolve().parent
TARGETS = {
    "eq3_search": "EQ3.cpp",
    "eq3_similarity": "eq3_similarity.cpp",
}

SETUP_TEMPLATE = '''\
import platform
import pybind11
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

if platform.system() == "Windows":
    extra_compile_args = ["/EHsc", "/bigobj", "/O2", "/std:c++17", "/utf-8"]
else:
    extra_compile_args = ["-O3", "-std=c++17"]

setup(
    name="{module}",
    ext_modules=[
        Pybind11Extension(
            "{module}",
            ["{source}"],
            include_dirs=[pybind11.get_include()],
            extra_compile_args=extra_compile_args,
            cxx_std=17,
        )
    ],
    cmdclass={{"build_ext": build_ext}},
    zip_safe=False,
)
'''


def loadable_artifacts(module):
    names = {module + suffix for suffix in importlib.machinery.EXTENSION_SUFFIXES}
    return sorted(p for p in CODE_DIR.iterdir() if p.name in names)


def can_import(module):
    probe = f"import sys; sys.path.insert(0, {str(CODE_DIR)!r}); import {module}"
    result = subprocess.run(
        [sys.executable, "-c", probe],
        capture_output=True,
        text=True,
    )
    if result.returncode == 0:
        return True, ""
    last_line = (result.stderr.strip().splitlines() or [""])[-1]
    return False, last_line


def ensure_pybind11():
    try:
        import pybind11  # noqa: F401
        return True
    except ImportError:
        pass

    print("  pybind11 未安装，正在安装")
    result = subprocess.run(
        [sys.executable, "-m", "pip", "install", "pybind11"],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print("✗ pybind11 安装失败")
        print(result.stderr)
        return False

    print("  ✓ pybind11 就绪")
    return True


def normalize_module_name(source_text, module):
    match = re.search(r"PYBIND11_MODULE\s*\(\s*([A-Za-z_]\w*)\s*,", source_text)
    if not match:
        return None, "源文件里找不到 PYBIND11_MODULE 定义"

    declared = match.group(1)
    if declared == module:
        return source_text, ""

    print(f"  源文件声明的模块名是 {declared}，改写为 {module}")
    patched = (
        source_text[:match.start(1)] + module + source_text[match.end(1):]
    )
    return patched, ""


def build(module, source_path):
    print(f"  源文件 {source_path}")

    source_text = source_path.read_text(encoding="utf-8")
    source_text, error = normalize_module_name(source_text, module)
    if source_text is None:
        print(f"✗ {error}")
        return False

    with tempfile.TemporaryDirectory(prefix=f"build_{module}_") as temp_dir:
        temp = Path(temp_dir)
        (temp / source_path.name).write_text(source_text, encoding="utf-8")
        (temp / "setup.py").write_text(
            SETUP_TEMPLATE.format(module=module, source=source_path.name),
            encoding="utf-8",
        )

        print("  开始编译")
        result = subprocess.run(
            [sys.executable, "setup.py", "build_ext", "--inplace"],
            cwd=temp,  # 不用 os.chdir，避免异常时污染进程状态
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
        )

        if result.returncode != 0:
            print("✗ 编译失败")
            print(result.stdout)
            print(result.stderr)
            return False

        suffixes = tuple(importlib.machinery.EXTENSION_SUFFIXES)
        built = [p for p in temp.iterdir() if p.name.startswith(module) and p.name.endswith(suffixes)]

        if not built:
            print("✗ 编译进程返回成功，但没有生成扩展文件")
            print(result.stdout)
            return False

        for artifact in built:
            destination = CODE_DIR / artifact.name
            shutil.copy2(artifact, destination)
            print(f"  ✓ 产物 {destination}")

    return True


def process(module, force):
    source_name = TARGETS[module]
    print("=" * 60)
    print(f"[{module}]")

    existing = loadable_artifacts(module)

    if not force and existing:
        ok, error = can_import(module)
        if ok:
            print(f"  ✓ 已有可用产物，跳过：{existing[0].name}")
            return True
        print(f"  产物存在但导入失败，重新编译：{error}")

    source_path = CODE_DIR / source_name

    if not source_path.exists():
        ok, error = can_import(module)
        if ok:
            print(f"  ✓ 源文件 {source_name} 不存在，但模块已能导入，跳过")
            return True
        print(f"✗ 源文件不存在，且模块无法导入：{source_path}")
        print(f"  导入错误：{error}")
        print(f"  请把 {source_name} 放到 {CODE_DIR}")
        return False

    if not build(module, source_path):
        return False

    ok, error = can_import(module)
    if not ok:
        print(f"✗ 编译产物无法导入：{error}")
        return False

    print(f"  ✓ import {module} 成功")
    return True


def main():
    parser = argparse.ArgumentParser(description="编译 EQ3 的 C++ 扩展模块")
    parser.add_argument("modules", nargs="*", default=[],
                        help=f"要编译的模块，默认全部；可选 {', '.join(TARGETS)}")
    parser.add_argument("--force", action="store_true", help="忽略已有产物，强制重编")
    args = parser.parse_args()

    unknown = [module for module in args.modules if module not in TARGETS]
    if unknown:
        parser.error(f"未知模块 {', '.join(unknown)}；可选 {', '.join(TARGETS)}")

    modules = args.modules or list(TARGETS)

    print(f"解释器 {sys.executable}")
    print(f"目标目录 {CODE_DIR}")
    print(f"扩展后缀 {importlib.machinery.EXTENSION_SUFFIXES[0]}")

    if not ensure_pybind11():
        return 1

    failed = [module for module in modules if not process(module, args.force)]

    print("=" * 60)
    if failed:
        print(f"✗ 失败：{', '.join(failed)}")
        return 1

    print("🎉 全部就绪，现在可以跑 run.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
