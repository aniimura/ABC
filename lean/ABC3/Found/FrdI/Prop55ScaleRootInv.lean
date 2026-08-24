/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop55BiratPfFun

/-!
# [FrdI] Proposition 5.5, (ii) —— 根の一斉倍化 `Σ_k` は 3 量を保つ

原文 (FrdI p.104):
> (ii) There is a natural equivalence of categories [compatible with the functors to the

★★`Prop55BiratPfFun.lean` の測定で、残るのは
「`scaleRootEquiv`(`Σ_k : 𝒞^pf ≌ 𝒞^pf`、在庫)を **birat 化へ降ろす**」1 段だと分かった。
在庫の `psiBiratCor`(`Corollary 4.10`)がその形をしており、要る入力は

  `hfwd : ∀ f, coaPreProp f → coaPreProp (Σ_k.map f)`

である。`coaPreProp = IsCoAngular ∧ IsPreStep` なので、
★**`Σ_k` が `rootDeg` / `rootBase` / `rootDiv` を保つ**ことを言えばよい。

## ★★なぜ保たれるか

`Σ_k` は `rtRootIso`(= `rootIso`)の共役であり、
`rootIso` は **代表元の射そのものを変えず添字だけ押し出す**(`rootIso_hom_mk`)。
在庫の変化則は

| 量 | `rootIso` での変化 |
|---|---|
| `pfDeg` | ★不変(`pfDeg_rootIso`) |
| `pfBase` | `Base a` と `Base b` で共役(`pfBase_rootIso`) |
| `pfDiv` | `Φ.map (Base a)` で押し出し、分母は `degFr a` 倍(`pfDiv_rootIso`) |

であり、`rootBase` の定義がちょうどその共役を打ち消す。
★要になるのは `rtLift_ext`(`rtExt A d ≫ rtLift A h = rtExt A t`)である。

## ★★★★★測って分かったこと(2026-08-25)—— `rootDiv` は**不変ではない**

`rootDeg` と `rootBase` は保たれる(本ファイルで閉じた)が、
★**`rootDiv` は `1/k` 倍される**:

```
rootDiv (Σ_k f) = Pf.divBy k (rootDiv f)
```

計算: `rootDiv` は `X.root * Y.root` で割る定義なので `Σ_k` では `k²·n·m` で割る。
一方 `pfDiv` は `pfDiv_rootIso` により押し出しで `k` 倍になるだけなので、
差し引き **`k` 倍だけ余分に割られる**。

★★**これは正しい** —— `⟨A,n⟩` は「`A` の `n` 乗根」であり
`Σ_k ⟨A,n⟩ = ⟨A, k·n⟩` は「`k·n` 乗根」なので、零因子が `1/k` になるのが当然である。
★したがって **`Σ_k` は「零因子を保つ関手」ではない**。
`Prop55BiratPf.lean` で圏同値として作れているのはそのためであって、
Frobenioid の射としてではない。

## ★★残る段取り(次に着手する人へ)

`psiBiratCor`(`Corollary 4.10`)の入力は
`hfwd : ∀ f, coaPreProp f → coaPreProp (Σ_k.map f)` であり、
`coaPreProp = IsCoAngular ∧ IsPreStep`。

* `IsPreStep = IsLinear ∧ IsBaseIsomorphism` …… `rootDeg` / `rootBase` の不変性で**済**
* `IsCoAngular` …… **純粋に圏論的な条件**(分解 `φ = γ ≫ β ≫ α` と `IsIso β`)であり、
  `Σ_k` が圏同値であることと `IsLinear` / `IsIsometric` / `IsPreStep` /
  `IsBaseIsomorphism` の保存から出る。
* `IsIsometric`(零因子が `0`)…… 上の**スケール則**から
  `rootDiv (Σ_k f) = 0 ↔ rootDiv f = 0`(`k • Pf.divBy k w = w`)で出る。

★つまり残るのは**スケール則 1 本とそれを使う 3 つの保存**である。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section ScaleRootInv

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}

/-- ★★★★**`Σ_k` は Frobenius 次数を保つ**。 -/
theorem rootDeg_scaleRootHom (k : ℕ+) {X Y : PfRootObj P F} (f : HomRoot P F X Y) :
    rootDeg (scaleRootHom (F := F) k f) = rootDeg f :=
  pfDeg_rtRootIso_inv X.obj Y.obj
    (show k * Y.root = k * Y.root from rfl) (show k * X.root = k * X.root from rfl) f

/-- ★★★★★**`Σ_k` は底を保つ**。

★`rootBase` は `rtExt` で共役を取る形なので、`rtLift_ext` で
`rtExt A (k·d) = rtExt A d ≫ rtLift` と割ると、`pfBase_rootIso` の共役がちょうど消える。 -/
theorem rootBase_scaleRootHom (k : ℕ+) {X Y : PfRootObj P F} (f : HomRoot P F X Y) :
    rootBase (scaleRootHom (F := F) k f) = rootBase f := by
  haveI iX : IsIso (P.Base (rtExt P F Y.obj X.root)) := (rtExt_frobType P F Y.obj X.root).2
  haveI iXk : IsIso (P.Base (rtExt P F Y.obj (k * X.root))) :=
    (rtExt_frobType P F Y.obj (k * X.root)).2
  set hA : k * Y.root = k * Y.root := rfl with hAdef
  set hB : k * X.root = k * X.root := rfl with hBdef
  set a := rtLift P F X.obj hA with hadef
  set b := rtLift P F Y.obj hB with hbdef
  haveI ib : IsIso (P.Base b) := (rtLift_frobType P F Y.obj hB).2
  -- ★`pfBase` の共役則(逆向き)
  have hcon : pfBase f ≫ P.Base b
      = P.Base a ≫ pfBase (scaleRootHom (F := F) k f) := by
    have h := pfBase_rootIso a (rtLift_frobType P F X.obj hA) b (rtLift_frobType P F Y.obj hB)
      (by rw [rtLift_degFr, rtLift_degFr]) (scaleRootHom (F := F) k f)
    rwa [show (rootIso (F := F) a (rtLift_frobType P F X.obj hA) b
        (rtLift_frobType P F Y.obj hB) (by rw [rtLift_degFr, rtLift_degFr])).hom
        (scaleRootHom (F := F) k f)
        = f from Iso.inv_hom_id_apply (rtRootIso P F X.obj Y.obj hA hB) f] at h
  -- ★`rtExt` を割る
  have hextX : P.Base (rtExt P F X.obj (k * Y.root))
      = P.Base (rtExt P F X.obj Y.root) ≫ P.Base a := by
    rw [← P.Base_comp, hadef, rtLift_ext]
  have hextY : P.Base (rtExt P F Y.obj (k * X.root))
      = P.Base (rtExt P F Y.obj X.root) ≫ P.Base b := by
    rw [← P.Base_comp, hbdef, rtLift_ext]
  show P.Base (rtExt P F X.obj (k * Y.root)) ≫ pfBase (scaleRootHom (F := F) k f)
      ≫ @inv _ _ _ _ (P.Base (rtExt P F Y.obj (k * X.root))) iXk
    = P.Base (rtExt P F X.obj Y.root) ≫ pfBase f
      ≫ @inv _ _ _ _ (P.Base (rtExt P F Y.obj X.root)) iX
  -- ★両辺に `Base (rtExt Y (k·X.root))` を掛けて比べる
  refine (cancel_mono (P.Base (rtExt P F Y.obj (k * X.root)))).mp ?_
  have hinvK : @inv _ _ _ _ (P.Base (rtExt P F Y.obj (k * X.root))) iXk
      ≫ P.Base (rtExt P F Y.obj (k * X.root)) = 𝟙 _ :=
    @IsIso.inv_hom_id _ _ _ _ (P.Base (rtExt P F Y.obj (k * X.root))) iXk
  have hinvX : @inv _ _ _ _ (P.Base (rtExt P F Y.obj X.root)) iX
      ≫ P.Base (rtExt P F Y.obj X.root) = 𝟙 _ :=
    @IsIso.inv_hom_id _ _ _ _ (P.Base (rtExt P F Y.obj X.root)) iX
  have hL : (P.Base (rtExt P F X.obj (k * Y.root)) ≫ pfBase (scaleRootHom (F := F) k f)
        ≫ @inv _ _ _ _ (P.Base (rtExt P F Y.obj (k * X.root))) iXk)
      ≫ P.Base (rtExt P F Y.obj (k * X.root))
      = P.Base (rtExt P F X.obj (k * Y.root)) ≫ pfBase (scaleRootHom (F := F) k f) :=
    (Category.assoc _ _ _).trans
      (congrArg (fun t => P.Base (rtExt P F X.obj (k * Y.root)) ≫ t)
        ((Category.assoc _ _ _).trans
          ((congrArg (fun t => pfBase (scaleRootHom (F := F) k f) ≫ t) hinvK).trans
            (Category.comp_id _))))
  have hR : (P.Base (rtExt P F X.obj Y.root) ≫ pfBase f
        ≫ @inv _ _ _ _ (P.Base (rtExt P F Y.obj X.root)) iX)
      ≫ P.Base (rtExt P F Y.obj (k * X.root))
      = P.Base (rtExt P F X.obj Y.root) ≫ pfBase f ≫ P.Base b := by
    have s : @inv _ _ _ _ (P.Base (rtExt P F Y.obj X.root)) iX
          ≫ (P.Base (rtExt P F Y.obj X.root) ≫ P.Base b)
        = (@inv _ _ _ _ (P.Base (rtExt P F Y.obj X.root)) iX
            ≫ P.Base (rtExt P F Y.obj X.root)) ≫ P.Base b := (Category.assoc _ _ _).symm
    rw [hextY, Category.assoc, Category.assoc, s, hinvX, Category.id_comp]
  rw [hL, hR, hextX, Category.assoc, ← hcon]
  rfl

end ScaleRootInv

/-! ### ★出典の紐付け -/

/-- ★★★★locator —— `Proposition 5.5, (ii)` の梱包に要る
「根の一斉倍化 `Σ_k` は `rootDeg` / `rootBase` を保つ」。 -/
def rootBase_scaleRootHom.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — Σ_k は Frobenius 次数と底を保つ",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
