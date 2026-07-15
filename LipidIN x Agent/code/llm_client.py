"""Pluggable multi-provider LLM client (Anthropic Claude / OpenAI-Codex / DeepSeek / custom).

可插拔的多提供方大模型客户端（Anthropic Claude / OpenAI-Codex / DeepSeek / 自定义）。

The user fills in the provider, API key and model. Anthropic uses the Messages API; everything else
uses the OpenAI-compatible /chat/completions schema (DeepSeek and most local gateways follow it).
Only the Python standard library is used (urllib), so there is no extra dependency to install.

用户填写提供方、API 密钥与模型。Anthropic 用 Messages API；其余走 OpenAI 兼容的 /chat/completions
（DeepSeek 及大多数本地网关都兼容）。仅用标准库（urllib），无需额外安装依赖。
"""

import os
import re
import json
import time
import urllib.request
import urllib.error


PROVIDER_DEFAULTS = {
    "anthropic": {"base_url": "https://api.anthropic.com", "model": "claude-sonnet-5",
                  "key_env": "ANTHROPIC_API_KEY", "kind": "anthropic"},
    "openai": {"base_url": "https://api.openai.com/v1", "model": "gpt-4o-mini",
               "key_env": "OPENAI_API_KEY", "kind": "openai"},
    "codex": {"base_url": "https://api.openai.com/v1", "model": "gpt-4o-mini",
              "key_env": "OPENAI_API_KEY", "kind": "openai"},
    "deepseek": {"base_url": "https://api.deepseek.com", "model": "deepseek-chat",
                 "key_env": "DEEPSEEK_API_KEY", "kind": "openai"},
    "custom": {"base_url": "", "model": "", "key_env": "LLM_API_KEY", "kind": "openai"},
}


class LLMClient:
    def __init__(self, provider="none", api_key=None, model=None, base_url=None,
                 max_tokens=1500, temperature=0.2, timeout=90):
        self.provider = (provider or "none").lower().strip()
        self.max_tokens = max_tokens
        self.temperature = temperature
        self.timeout = timeout

        if self.provider in ("none", ""):
            self.enabled = False
            self.kind = None
            self.model = None
            self.base_url = None
            self.api_key = None
            return

        defaults = PROVIDER_DEFAULTS.get(self.provider, PROVIDER_DEFAULTS["custom"])
        self.kind = defaults["kind"]
        self.model = model or defaults["model"]
        self.base_url = (base_url or defaults["base_url"]).rstrip("/")
        self.api_key = api_key or os.environ.get(defaults["key_env"]) or os.environ.get("LLM_API_KEY")
        self.enabled = bool(self.api_key) and bool(self.base_url) and bool(self.model)

    @classmethod
    def from_dict(cls, config):
        config = config or {}
        return cls(
            provider=config.get("provider", "none"),
            api_key=config.get("api_key"),
            model=config.get("model"),
            base_url=config.get("base_url"),
            max_tokens=config.get("max_tokens", 1500),
            temperature=config.get("temperature", 0.2),
            timeout=config.get("timeout", 90),
        )

    def is_enabled(self):
        return bool(self.enabled)

    # ------------------------------------------------------------------ #
    def chat(self, system, user, retries=2, json_mode=False, max_tokens=None):
        """Send one system+user turn and return the assistant text. Raises on failure.

        发送一次 system+user 对话并返回助手文本，失败时抛异常。
        """
        if not self.enabled:
            raise RuntimeError("LLM client is not configured / 大模型客户端未配置")

        last_error = None
        for attempt in range(retries + 1):
            try:
                if self.kind == "anthropic":
                    return self._chat_anthropic(system, user, json_mode, max_tokens)
                return self._chat_openai(system, user, json_mode, max_tokens)
            except Exception as e:
                last_error = e
                if attempt < retries:
                    time.sleep(1.5 * (attempt + 1))
        raise RuntimeError(f"LLM request failed after retries / 大模型请求多次失败: {last_error}")

    def chat_json(self, system, user, retries=2):
        """Ask for JSON and parse it robustly: JSON mode, a bigger token budget, corrective retries.

        请求并稳健解析 JSON：启用 JSON 模式、放大 token 预算、失败时带纠正提示重试。
        """
        max_tokens = max(self.max_tokens, 2048)
        instruction = ("\n\nReturn ONLY a single valid JSON object (no prose, no markdown, no code fences). "
                       "Start your reply with '{' and end with '}'.")
        attempts = max(2, retries + 1)
        last_text = ""
        last_error = None
        for i in range(attempts):
            extra = instruction
            if i > 0:
                extra += (" Your previous reply was empty or not valid JSON. "
                          "Output the JSON object only, nothing else.")
            try:
                text = self.chat(system, user + extra, retries=1, json_mode=True, max_tokens=max_tokens)
            except Exception as e:
                last_error = e
                continue
            last_text = text or ""
            if last_text.strip():
                try:
                    return extract_json(last_text)
                except ValueError as e:
                    last_error = e
        detail = repr(last_text[:200]) if last_text else (str(last_error) if last_error else "empty")
        raise ValueError(f"No JSON in LLM response after {attempts} attempts "
                         f"/ 多次尝试仍未取得有效 JSON: {detail}")

    # ------------------------------------------------------------------ #
    def _post(self, url, headers, payload):
        data = json.dumps(payload).encode("utf-8")
        request = urllib.request.Request(url, data=data, headers=headers, method="POST")
        try:
            with urllib.request.urlopen(request, timeout=self.timeout) as response:
                return json.loads(response.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            body = e.read().decode("utf-8", errors="replace")
            raise RuntimeError(f"HTTP {e.code}: {body[:400]}")

    def _chat_anthropic(self, system, user, json_mode=False, max_tokens=None):
        url = f"{self.base_url}/v1/messages"
        headers = {
            "x-api-key": self.api_key,
            "anthropic-version": "2023-06-01",
            "content-type": "application/json",
        }
        messages = [{"role": "user", "content": user}]
        if json_mode:
            # Prefill the assistant turn with "{" so the model is forced to emit a JSON object.
            # 用 "{" 预填助手轮，强制模型输出 JSON 对象。
            messages.append({"role": "assistant", "content": "{"})
        payload = {
            "model": self.model,
            "max_tokens": max_tokens or self.max_tokens,
            "temperature": self.temperature,
            "system": system,
            "messages": messages,
        }
        result = self._post(url, headers, payload)
        parts = result.get("content", [])
        texts = [p.get("text", "") for p in parts if isinstance(p, dict) and p.get("type") == "text"]
        text = "".join(texts).strip()
        if json_mode and text and not text.lstrip().startswith("{"):
            text = "{" + text  # restore the prefill that the API does not echo back
        return text

    def _chat_openai(self, system, user, json_mode=False, max_tokens=None):
        url = f"{self.base_url}/chat/completions"
        headers = {
            "Authorization": f"Bearer {self.api_key}",
            "content-type": "application/json",
        }
        payload = {
            "model": self.model,
            "max_tokens": max_tokens or self.max_tokens,
            "temperature": self.temperature,
            "messages": [
                {"role": "system", "content": system},
                {"role": "user", "content": user},
            ],
        }
        # Native JSON mode for providers known to support it (OpenAI, DeepSeek).
        # 对已知支持的提供方（OpenAI、DeepSeek）启用原生 JSON 模式。
        if json_mode and self.provider in ("openai", "codex", "deepseek"):
            payload["response_format"] = {"type": "json_object"}
        result = self._post(url, headers, payload)
        return _extract_openai_content(result)


def _extract_openai_content(result):
    """Pull assistant text from an OpenAI-compatible response, tolerating list content and reasoning-only replies.

    从 OpenAI 兼容响应中取出助手文本，兼容列表型 content 与只有 reasoning 的回复。
    """
    choices = result.get("choices") or []
    if not choices:
        # Some gateways surface errors as a 200 body with an 'error' field.
        # 某些网关把错误放在 200 响应的 'error' 字段里。
        if "error" in result:
            raise RuntimeError(f"LLM error: {str(result['error'])[:300]}")
        return ""
    message = choices[0].get("message") or {}
    content = message.get("content")
    if isinstance(content, list):
        content = "".join(
            part.get("text", "") for part in content
            if isinstance(part, dict) and part.get("type") in ("text", "output_text", None)
        )
    if not content:
        # Reasoning models may leave 'content' empty and put text in 'reasoning_content'.
        # 推理模型可能把正文留空，而把文本放进 'reasoning_content'。
        content = message.get("reasoning_content") or ""
    return (content or "").strip()


def extract_json(text):
    """Pull the first JSON object/array out of an LLM reply, tolerating ```json fences and prose.

    从大模型回复里取出第一个 JSON 对象/数组，容忍 ```json 围栏和多余文字。
    """
    if text is None or not str(text).strip():
        raise ValueError("empty LLM response / 大模型响应为空")
    fenced = re.search(r"```(?:json)?\s*(.*?)```", text, flags=re.DOTALL)
    candidate = fenced.group(1).strip() if fenced else text.strip()
    try:
        return json.loads(candidate)
    except json.JSONDecodeError:
        pass
    for opener, closer in (("{", "}"), ("[", "]")):
        start = candidate.find(opener)
        end = candidate.rfind(closer)
        if start != -1 and end != -1 and end > start:
            try:
                return json.loads(candidate[start:end + 1])
            except json.JSONDecodeError:
                continue
    raise ValueError(f"No JSON found in LLM response / 未在大模型响应中找到 JSON: {text[:200]}")
