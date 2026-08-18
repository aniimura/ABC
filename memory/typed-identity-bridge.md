---
name: typed-identity-bridge
description: 台が同じで綴りが違う型の橋は、instance を足すより型付き恒等関数を定義する方が安全
metadata:
  type: feedback
---

Lean で「台は同じだが綴りが違う」型(`(𝟙_ …).obj t` と `Γ(X,U)`、`ModuleCat` の係数環の `CommRingCat` 綴りと `RingCat` 綴り、など)の間を渡るときは、**型付き恒等関数**を定義する:

```lean
def unitVal (x : ((𝟙_ (PresheafModulesOn X U)).obj t : Type u)) : (Γ(X,U) : Type u) := x
```

すると `unitVal (a • b) = unitVal a * unitVal b` や `unitVal 1 = 1` が **`rfl` で通る**ことが多い。

**Why:** 2026-08-18、B2 の `unitEnd` で 5 通り試して 4 通りが失敗した。特に `inferInstanceAs` で `CommRing` instance を足す手は、**`SMul` の経路が 2 本になり `rw` が当たらなくなる**という新しい問題を作った。型付き恒等関数は既存の経路と競合しない。

**How to apply:** 型注釈 `(x : T)` は**推論された型を変えない**ので橋にならない(本 session で 3 回踏んだ)。`have h : T := x` か、再利用するなら `def` で恒等関数を置く。橋の補題(`_smul`, `_one` など)は**定義と同じファイル**に置く——別ファイルだと `show` で分解する際に型が決まらない。[[ring-instance-two-paths]] の逃げ道の改良版。
