---
name: ring-instance-two-paths
description: mathlib の前層加群で `Ring` インスタンスが 2 経路で付き、simp/rw が発火しなくなる罠と塞ぎ方
metadata:
  type: project
---

`↑((R ⋙ forget₂ CommRingCat RingCat).obj Z)` には `Ring` インスタンスが **2 経路**で付く。

| 経路 | 出どころ | 使う側 |
|---|---|---|
| (a) | `RingCat` の構造 | `PresheafOfModules.freeObj` / `ModuleCat.restrictScalars` |
| (b) | `CommRing.toRing`(`Mathlib/Algebra/Category/ModuleCat/Presheaf/Monoidal.lean` の `CommRing` インスタンス) | 前層加群のテンソル積 |

**Why:** (a) と (b) は defeq だが**構文が違う**。両方が混ざった項では `simp`/`rw` が
`CategoryTheory.comp_apply`(`@[simp]` 付き)にすら当たらない。
2026-08-18 に `free (yoneda V) ⊗ free (yoneda W) ≅ free (yoneda (V ⊓ W))` の
naturality を書こうとして 15 回踏んだ。

**How to apply:** 項の**構成**は通る(型検査は defeq を見る)ので、
詰まるのは必ず**証明**側である。塞ぎ方は
(1) `letI` で片方に固定、(2) `ModuleCat.of` で束ね直す、
(3) **`show` で両辺を明示的に書き下してから `rw`**——(3) が本命。
[[measure-before-deciding]] と同じく、まず「どちらの経路の項か」を実測すること。
これは「インスタンス束縛子は型の書き方の違いをまたげない」(探索の失敗)の**兄弟**で、
こちらは**書き換えの失敗**である。
