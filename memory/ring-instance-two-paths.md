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

## ★★★★★★2026-08-19 追加 —— **書き換えられないなら、書き換えずに主張する**

同じ罠が `Scheme.Modules` のファイバーで出た。`arcFiber p L` の `ℂ` 作用は
`restrictScalars (ΓSpecIso).inv` を通しており、係数環 `Γ(Spec ℂ, ⊤)` の作用と
**定義的に等しいが綴りが違う**。

★試して**全部落ちた**もの:

| 手 | 結果 |
|---|---|
| `rw [ModuleCat.restrictScalars.smul_def]` | ★パターン不一致 |
| `show ... (c' • v) ...` | ★`HSMul Γ ↥(arcFiber p L)` のインスタンスが無い |
| 型上書き `(v : ↥(arcFiber p L))` | ★★**効かない**——`v` の推論型は変わらない |
| `show T from v` | ★`have this := v; this` になり型が変わらない |

★★★通ったのはこれだけ:

    have h : topMap e (c • v) = ΓSpecIso.inv.hom c * topMap e v :=
      topMap_smul e (ΓSpecIso.inv.hom c) v

★`have` の**型はゴール側の綴り**で書き、証明項は**別の綴り**の補題を渡す。

**Why:** ★★★★★**インスタンス探索は構文的な型を見るが、項の型検査は定義的等しさを見る。**
`rw` / `simp` / `show` はすべて構文照合を経由するので落ちるが、
`have ... := <別綴りの補題>` は `isDefEq` だけを通るので受理される。

**How to apply:**
- 2 経路の綴りで詰まったら、**まず `have` でゴール側の型を書き、別綴りの補題を `exact` で渡す**。
- `rw [map_zero]` 等が `presheaf` の二重路(`TopCat.Presheaf` vs `_ᵒᵖ ⥤ _`)で落ちたら、
  `(congrArg (fun y => f y) h).trans (map_zero _)` に置き換える。
- 関連: [[lean-rebind-morphisms-clean-types]] / [[typed-identity-bridge]] /
  [[defer-fixing-coordinates]]
