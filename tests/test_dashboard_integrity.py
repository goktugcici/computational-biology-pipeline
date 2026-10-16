import re
import os
import sys

def test_html_integrity():
    html_path = "dashboard/index.html"
    assert os.path.exists(html_path), "dashboard/index.html not found!"
    with open(html_path, "r", encoding="utf-8") as f:
        html = f.read()

    print(f"[1] HTML Size: {len(html)} characters")

    # 1. Script blocks
    script_blocks = re.findall(r"<script(?:\s+[^>]*)?>(.*?)</script>", html, re.DOTALL)
    print(f"[2] Script blocks found: {len(script_blocks)}")
    full_js = "\n".join(script_blocks)

    # 2. Extract DOM IDs
    dom_ids = set(re.findall(r'id=["\']([^"\']+)["\']', html))
    print(f"[3] Unique DOM IDs defined in HTML: {len(dom_ids)}")

    # 3. Check getElementById calls
    get_elem_ids = set(re.findall(r'getElementById\(["\']([^"\']+)["\']\)', full_js))
    print(f"[4] Unique getElementById calls in JS: {len(get_elem_ids)}")

    missing_ids = [gid for gid in get_elem_ids if gid not in dom_ids]
    if missing_ids:
        print(f"  [!] WARNING: getElementById targets not found in DOM: {missing_ids}")
    else:
        print("  [OK] All getElementById targets exist in HTML DOM!")

    # 4. Check all tabs and tab switching
    tab_targets = set(re.findall(r'switchTab\(["\']([^"\']+)["\']\)', html))
    print(f"[5] switchTab targets ({len(tab_targets)}): {sorted(list(tab_targets))}")

    missing_tabs = []
    for tab in tab_targets:
        tab_section = f"tab-{tab}"
        tab_btn = f"tabBtn-{tab}"
        if tab_section not in dom_ids:
            missing_tabs.append(f"Missing section ID: {tab_section}")
        if tab_btn not in dom_ids:
            missing_tabs.append(f"Missing button ID: {tab_btn}")

    if missing_tabs:
        print(f"  [!] ERROR: Missing tabs: {missing_tabs}")
    else:
        print(f"  [OK] All {len(tab_targets)} switchTab targets have matching tab-* and tabBtn-* elements in DOM!")

    # 5. Check inline handlers
    handlers = re.findall(r'''(on[a-z]+)\s*=\s*["']([^"']+)["']''', html)
    print(f"[6] Inline Event Handlers ({len(handlers)})")

    # Extract all functions declared in JS
    declared_funcs = set(re.findall(r'function\s+([a-zA-Z0-9_]+)\s*\(', full_js))
    declared_consts = set(re.findall(r'(?:const|let|var)\s+([a-zA-Z0-9_]+)\s*=\s*(?:function|\([^)]*\)\s*=>)', full_js))
    all_defined = declared_funcs.union(declared_consts)

    builtins = {'alert', 'confirm', 'prompt', 'console', 'window', 'document', 'open', 'print', 'focus', 'blur'}
    keywords = {'if', 'for', 'while', 'switch', 'return', 'let', 'const', 'var', 'event'}

    missing_handler_calls = []
    for ev, code in handlers:
        calls = re.findall(r"([a-zA-Z0-9_]+)\s*\(", code)
        for c in calls:
            if c not in all_defined and c not in builtins and c not in keywords:
                missing_handler_calls.append((ev, code, c))

    if missing_handler_calls:
        print(f"  [!] ERROR: Handlers calling undeclared functions: {missing_handler_calls}")
    else:
        print("  [OK] All inline handlers call properly defined JavaScript functions!")

    # 6. Verify external libraries
    scripts = re.findall(r'<script\s+src=["\']([^"\']+)["\']', html)
    print(f"[7] External dependencies: {scripts}")

    is_all_clean = len(missing_ids) == 0 and len(missing_tabs) == 0 and len(missing_handler_calls) == 0
    print("\n[RESULT]:", "PASSED" if is_all_clean else "FAILED")
    return is_all_clean

if __name__ == "__main__":
    success = test_html_integrity()
    sys.exit(0 if success else 1)
