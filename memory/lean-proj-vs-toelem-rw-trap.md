---
name: lean-proj-vs-toelem-rw-trap
description: `rw`/`simp` が「型が正しくない」と言って止まる典型——defeq だが構文が違う 2 つの綴りが混ざった目標
metadata:
  type: feedback
---

`rw` / `simp` が **「Did not find an occurrence」＋「The target expression is not type-correct under the `instances` transparency level」** と言って止まったら、
★**defeq だが構文の違う 2 つの綴りが目標に混ざっていないか**を最初に疑う。

**Why:** `rw` は書き換え後の項を `instances` 透明度で型検査する。`def` を展開しないと一致しない 2 つの綴り（例: `P.proj.obj X` と `(P.toElem.obj X).base`、`P.proj` は `P.toElem ⋙ ElemFrobCat.proj` という `def`）が混ざっていると、書き換え後の項がその透明度で型付かず、`rw` は失敗する。**インスタンス引数の違いでも `rw` の適用位置でもない。**

**How to apply:**
1. エラー本文の「Full error: Application type mismatch」を読む。**そこに 2 つの綴りが並んで出る**——それが原因。
2. 対処は **`rw` を使わないこと**。`calc` と `congrArg` / `Category.assoc` を**項として**書けば、既定の透明度で defeq が通る。
3. 目標に `inv f` などインスタンスを担ぐ項があると、型が確定するまで TC が走らず `failed to synthesize` になる。`obtain ⟨a, ha⟩ : ∃ a : 望む型, a = 元の項 := ⟨_, rfl⟩` で**綴りの決まった変数を先に導入**すると通る（`Prop16.lean` の分類表 #3 と同じ手）。

★2026-08-20 実測: [FrdI] `Proposition 1.6, (ii)` の `plBkEquiv` はこれで 1 セッション止まっていた（前任の候補 A「`inv` のインスタンス」・候補 B「`rw` の適用位置」はどちらも外れ）。原因を特定したら 30 分で閉じた。

関連: [[frdi-s5-s6-blockers]]
