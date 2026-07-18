# -*- coding: utf-8 -*-
"""
兼容性垫片(shim)。

EQ2_batchsearchv2.compile_cpp_module() 会先 `import Cpp_batch_search` 作为
“C++ 搜库模块是否就绪”的门槛检查；只有该 import 成功(返回 True)，EQ3_batch 才会
继续。真正执行搜库的是预编译好的 eq3_search(.pyd)，与本模块名无关。

本机没有 C++ 编译器(cl.exe/g++ 均缺失)，无法用 build_cpp.py 现场编译出
Cpp_batch_search.pyd；但 eq3_search / eq3_similarity 的 .pyd 已存在且可正常 import。
因此这里提供一个纯 Python 垫片：转导出 eq3_search，让门槛检查通过即可。
"""
from eq3_search import *  # noqa: F401,F403

try:
    import eq3_search as _eq3_search  # 保留原模块引用，便于需要时访问
except Exception:  # pragma: no cover
    _eq3_search = None
