---
name: namespace-shadows-mathlib
description: プロジェクトの namespace が mathlib の同名 namespace を隠す——ABC3.Found.NumberField があると open NumberField が mathlib を指さない
metadata:
  type: feedback
---

`ABC3/Found/NumberField/PrimeDivisorsOfValues.lean` は `namespace ABC3.Found.NumberField` を
宣言している。そのため **`namespace ABC3.Found.NF` の内側で `open NumberField` と書くと、
Lean は相対解決で `ABC3.Found.NumberField` を先に見つけ、mathlib の `NumberField` が隠れる**。
症状は「`𝓞` が unknown identifier」で、autoImplicit が `𝓞` を型変数に束縛するため
エラーが**遠くの行に化けて出る**。

**Why:** 2026-08-21、`Found/NumberField/SplitInfinite.lean` でこれに 3 往復とられた。
同じディレクトリの `SplitCount.lean` は動いていたので「import が足りない」と誤診しかけた
——**動いていたのは、その時点で `PrimeDivisorsOfValues` を import していなかったから**である。

**How to apply:**
* `ABC3/Found/NumberField/` 配下では `open _root_.NumberField` / `open scoped _root_.NumberField` と書く。
* mathlib の定数を参照するときも `_root_.NumberField.Ideal.…` と書く。
* ★症状が「見当違いの行の unknown identifier」なら、まず **最小再現**（import と open だけの
  probe ファイル）を作って二分すること。→ [[python-w-truncates]]（probe が空でないことも確かめる）
