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

/-! ## ★2. `Pf.divBy` の道具 -/

/-- ★2 段の割り算は 1 段。 -/
theorem Pf.divBy_divBy {M : Type w} [AddCommMonoid M] (k l : ℕ+) (x : Pf M) :
    Pf.divBy k (Pf.divBy l x) = Pf.divBy (k * l) x := by
  induction x using Pf.inductionOn with | _ m a =>
  rw [Pf.divBy_mk, Pf.divBy_mk, Pf.divBy_mk, mul_assoc]

/-- ★`Pf.map` は割り算と可換。 -/
theorem Pf.map_divBy {M N : Type w} [AddCommMonoid M] [AddCommMonoid N] (f : M →+ N)
    (k : ℕ+) (x : Pf M) : Pf.map f (Pf.divBy k x) = Pf.divBy k (Pf.map f x) := by
  induction x using Pf.inductionOn with | _ m a =>
  rw [Pf.divBy_mk, Pf.map_mk, Pf.map_mk, Pf.divBy_mk]

/-! ## ★3. `rootDiv` の**スケール則** -/

/-- ★★★★★★**`Σ_k` は零因子を `1/k` 倍する**(不変ではない)。

★`rootDiv` は `X.root * Y.root` で割る定義なので `Σ_k` では `k²·n·m` で割るが、
`pfDiv` は押し出しで `k` 倍になるだけなので、差し引き `k` 倍だけ余分に割られる。
★★`⟨A,n⟩` が「`A` の `n` 乗根」であることを思えばこれが正しい ——
`Σ_k ⟨A,n⟩ = ⟨A, k·n⟩` は「`k·n` 乗根」だから零因子は `1/k` になる。 -/
theorem rootDiv_scaleRootHom (k : ℕ+) {X Y : PfRootObj P F} (f : HomRoot P F X Y) :
    rootDiv (scaleRootHom (F := F) k f) = Pf.divBy k (rootDiv f) := by
  set hA : k * Y.root = k * Y.root := rfl with hAdef
  set hB : k * X.root = k * X.root := rfl with hBdef
  set a := rtLift P F X.obj hA with hadef
  set b := rtLift P F Y.obj hB with hbdef
  set w := Pf.map (Φ.map (P.Base a)) (pfDiv (scaleRootHom (F := F) k f)) with hwdef
  -- ★`pfDiv` の変化則(`k` 倍になる)
  have hdiv : pfDiv f = Pf.divBy k w := by
    have h := pfDiv_rootIso a (rtLift_frobType P F X.obj hA) b (rtLift_frobType P F Y.obj hB)
      (by rw [rtLift_degFr, rtLift_degFr]) (scaleRootHom (F := F) k f)
    rw [show (rootIso (F := F) a (rtLift_frobType P F X.obj hA) b
        (rtLift_frobType P F Y.obj hB) (by rw [rtLift_degFr, rtLift_degFr])).hom
        (scaleRootHom (F := F) k f)
        = f from Iso.inv_hom_id_apply (rtRootIso P F X.obj Y.obj hA hB) f,
      show P.degFr a = k from rtLift_degFr P F X.obj hA] at h
    exact h
  -- ★押し出しを 2 段に割る
  have hmap : ∀ u : Pf (Φ.val (P.toElem.obj (rtObj P F X.obj (k * Y.root))).base),
      Pf.map (Φ.map (P.Base (rtExt P F X.obj (k * Y.root)))) u
        = Pf.map (Φ.map (P.Base (rtExt P F X.obj Y.root))) (Pf.map (Φ.map (P.Base a)) u) := by
    intro u
    have hbase : P.Base (rtExt P F X.obj (k * Y.root))
        = P.Base (rtExt P F X.obj Y.root) ≫ P.Base a := by
      rw [← P.Base_comp, hadef, rtLift_ext]
    induction u using Pf.inductionOn with | _ m s =>
    rw [Pf.map_mk, Pf.map_mk, Pf.map_mk, hbase, MonoidOn.map_comp]
  show Pf.divBy ((k * X.root) * (k * Y.root))
      (Pf.map (Φ.map (P.Base (rtExt P F X.obj (k * Y.root))))
        (pfDiv (scaleRootHom (F := F) k f)))
    = Pf.divBy k (Pf.divBy (X.root * Y.root)
        (Pf.map (Φ.map (P.Base (rtExt P F X.obj Y.root))) (pfDiv f)))
  rw [hmap, ← hwdef, hdiv, Pf.map_divBy, Pf.divBy_divBy, Pf.divBy_divBy]
  exact congrArg (fun t : ℕ+ => Pf.divBy t
      (Pf.map (Φ.map (P.Base (rtExt P F X.obj Y.root))) w))
    (by simp only [mul_comm, mul_assoc, mul_left_comm] :
      (k * X.root) * (k * Y.root) = k * (X.root * Y.root) * k)

/-! ## ★4. `Σ_k` は `Definition 1.2` の型を保つ(`IsCoAngular` を除く) -/

/-- ★`Σ_k` は linear を保つ。 -/
theorem isLinear_scaleRootHom (k : ℕ+) {X Y : PfRootObj P F} {f : HomRoot P F X Y}
    (h : IsLinear (pfRootPre P F) f) :
    IsLinear (pfRootPre P F) (scaleRootHom (F := F) k f) := by
  show rootDeg (scaleRootHom (F := F) k f) = 1
  rw [rootDeg_scaleRootHom]
  exact h

/-- ★`Σ_k` は base-isomorphism を保つ。 -/
theorem isBaseIsomorphism_scaleRootHom (k : ℕ+) {X Y : PfRootObj P F} {f : HomRoot P F X Y}
    (h : IsBaseIsomorphism (pfRootPre P F) f) :
    IsBaseIsomorphism (pfRootPre P F) (scaleRootHom (F := F) k f) := by
  show IsIso (rootBase (scaleRootHom (F := F) k f))
  rw [rootBase_scaleRootHom]
  exact h

/-- ★`Σ_k` は pre-step を保つ。 -/
theorem isPreStep_scaleRootHom (k : ℕ+) {X Y : PfRootObj P F} {f : HomRoot P F X Y}
    (h : IsPreStep (pfRootPre P F) f) :
    IsPreStep (pfRootPre P F) (scaleRootHom (F := F) k f) :=
  ⟨isLinear_scaleRootHom k h.1, isBaseIsomorphism_scaleRootHom k h.2⟩

/-- ★★**`Σ_k` は isometric を保つ**(零因子は `1/k` 倍されるが `0` は `0`)。 -/
theorem isIsometric_scaleRootHom (k : ℕ+) {X Y : PfRootObj P F} {f : HomRoot P F X Y}
    (h : IsIsometric (pfRootPre P F) f) :
    IsIsometric (pfRootPre P F) (scaleRootHom (F := F) k f) := by
  show rootDiv (scaleRootHom (F := F) k f) = 0
  rw [rootDiv_scaleRootHom, show rootDiv f = 0 from h, Pf.divBy_zero]
  rfl

/-- ★★**逆も成り立つ** —— `k • Pf.divBy k w = w` だから。 -/
theorem isIsometric_of_scaleRootHom (k : ℕ+) {X Y : PfRootObj P F} {f : HomRoot P F X Y}
    (h : IsIsometric (pfRootPre P F) (scaleRootHom (F := F) k f)) :
    IsIsometric (pfRootPre P F) f := by
  have h2 : Pf.divBy k (rootDiv f) = 0 := by
    rw [← rootDiv_scaleRootHom]
    exact h
  have h3 := congrArg (fun t : Pf (Φ.val (P.toElem.obj X.obj).base) => ((k : ℕ+) : ℕ) • t) h2
  rw [Pf.nsmul_divBy, smul_zero] at h3
  exact h3

/-- ★`Σ_k` は linear を反射する。 -/
theorem isLinear_of_scaleRootHom (k : ℕ+) {X Y : PfRootObj P F} {f : HomRoot P F X Y}
    (h : IsLinear (pfRootPre P F) (scaleRootHom (F := F) k f)) :
    IsLinear (pfRootPre P F) f := by
  have : rootDeg (scaleRootHom (F := F) k f) = 1 := h
  rwa [rootDeg_scaleRootHom] at this

/-- ★`Σ_k` は base-isomorphism を反射する。 -/
theorem isBaseIsomorphism_of_scaleRootHom (k : ℕ+) {X Y : PfRootObj P F}
    {f : HomRoot P F X Y}
    (h : IsBaseIsomorphism (pfRootPre P F) (scaleRootHom (F := F) k f)) :
    IsBaseIsomorphism (pfRootPre P F) f := by
  have : IsIso (rootBase (scaleRootHom (F := F) k f)) := h
  rwa [rootBase_scaleRootHom] at this

/-- ★`Σ_k` は pre-step を反射する。 -/
theorem isPreStep_of_scaleRootHom (k : ℕ+) {X Y : PfRootObj P F} {f : HomRoot P F X Y}
    (h : IsPreStep (pfRootPre P F) (scaleRootHom (F := F) k f)) :
    IsPreStep (pfRootPre P F) f :=
  ⟨isLinear_of_scaleRootHom k h.1, isBaseIsomorphism_of_scaleRootHom k h.2⟩

end ScaleRootInv

/-! ### ★出典の紐付け -/

/-- ★★★★locator —— `Proposition 5.5, (ii)` の梱包に要る
「根の一斉倍化 `Σ_k` は `rootDeg` / `rootBase` を保つ」。 -/
def rootBase_scaleRootHom.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — Σ_k は Frobenius 次数と底を保つ",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
