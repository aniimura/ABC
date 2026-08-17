import ABC3.Found.FrdI.Prop44Pre

/-!
# [FrdI] Proposition 4.4 の図式の中段 —— `𝒞^birat → 𝔽_{Φ^gp}` の零因子

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.83。

原文 (FrdI p.83):
> where the functors between elementary Frobenioids are those induced by the natural

★`Prop44.lean` は図式の**下段** `𝒞^birat → 𝔽_{0_𝒟}` までを作っている。
★★**中段 `𝒞^birat → 𝔽_{Φ^gp}` が無い**ので、`Proposition 4.4, (iii)` の
`Φ^birat`(その像として定まる)も書けない —— ★**ここが §4 の律速である**。

## ★★中段の零因子の式(紙の上で確定、2026-08-17)

★代表元 `(a : A′ → A, φ : A′ → B)` は形式的な合成 `a⁻¹ ≫ φ` を表す。
`𝔽_{Φ^gp}` の合成則(`Div_comp`)から

```
Div(a⁻¹)      = − Φ^gp(Base a)⁻¹ (Div a)
Div(a⁻¹ ≫ φ) = Φ^gp(Base a)⁻¹ (Div φ) − degFr(φ) · Φ^gp(Base a)⁻¹ (Div a)
```

★★**`degFr φ` の係数を落とさないこと** —— `Div_comp` の第 2 項が
`degFr φ • Div ψ` だからである(初稿でここを落としかけた)。

## ★★代表元の取り替えで不変であること

★遷移 `c`(co-angular pre-step、`c ≫ a_Z = a_W`)で `φ ↦ c ≫ φ` としたとき、
`Div c` の項が **2 か所に現れて打ち消す** ——
第 1 項の `degFr φ • Div c` と、第 2 項の `degFr φ · (Div a_W の Div c 成分)` である。
★★これが「不変量である」ことの中身で、**`a` が pre-step(次数 1)**が効く。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

/-! ## ★0. `toGp` は加法的 -/

variable {M : Type w} [AddCommMonoid M]

@[simp] theorem toGp_zero : toGp M 0 = 0 := by
  show (AddLocalization.addMonoidOf (⊤ : AddSubmonoid M)) 0 = 0
  exact map_zero _

theorem toGp_add (x y : M) : toGp M (x + y) = toGp M x + toGp M y := by
  show (AddLocalization.addMonoidOf (⊤ : AddSubmonoid M)) (x + y)
    = (AddLocalization.addMonoidOf (⊤ : AddSubmonoid M)) x
      + (AddLocalization.addMonoidOf (⊤ : AddSubmonoid M)) y
  exact map_add _ x y

theorem toGp_nsmul (n : ℕ) (x : M) : toGp M (n • x) = n • toGp M x := by
  induction n with
  | zero => simp
  | succ k ih => rw [succ_nsmul, toGp_add, ih, succ_nsmul]

/-! ## ★1. 代表元が定める `Φ^gp(A_𝒟)` の元 -/

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

variable {P} in
/-- ★★**代表元 `(a, φ)` が定める `Φ^gp(A_𝒟)` の元** —— 形式的な合成 `a⁻¹ ≫ φ` の零因子。

★`sliceBaseOf` と同じく **`IsIso` を仮引数に割る**(`Prop44.lean` の測定)。 -/
noncomputable def sliceDivGpOf {A B A' : C} (a : A' ⟶ A) (_ha : IsIso (P.Base a))
    (φ : A' ⟶ B) : Gp (Φ.val (P.toElem.obj A).base) :=
  haveI := _ha
  gpMap _ (Φ.map (inv (P.Base a)))
    (toGp _ (P.Div φ) - ((P.degFr φ : ℕ+) : ℕ) • toGp _ (P.Div a))

variable {P} in
theorem sliceDivGpOf_eq {A B A' : C} (a : A' ⟶ A) (ha : IsIso (P.Base a)) (φ : A' ⟶ B) :
    sliceDivGpOf a ha φ = haveI := ha
      gpMap _ (Φ.map (inv (P.Base a)))
        (toGp _ (P.Div φ) - ((P.degFr φ : ℕ+) : ℕ) • toGp _ (P.Div a)) := rfl

variable {P} in
/-- ★`IsIso` は `Prop` なので、射が等しければ値も等しい。 -/
theorem sliceDivGpOf_congr {A B A' : C} {a a' : A' ⟶ A} (h : a = a') (ha : IsIso (P.Base a))
    (ha' : IsIso (P.Base a')) (φ : A' ⟶ B) :
    sliceDivGpOf a ha φ = sliceDivGpOf a' ha' φ := by
  subst h; rfl

/-! ## ★2. ★★★代表元の取り替えで不変 -/

variable (Φ) in
/-- ★`gpMap ∘ gpMap` を `Φ.map` の合成にまとめる。 -/
theorem gpMap_phi_comp {A B E : D} (β : A ⟶ B) (α : B ⟶ E) (x : Gp (Φ.val E)) :
    gpMap _ (Φ.map β) (gpMap _ (Φ.map α) x) = gpMap _ (Φ.map (β ≫ α)) x := by
  have h : (Φ.map β).comp (Φ.map α) = Φ.map (β ≫ α) := by
    ext y
    exact (Φ.map_comp α β y).symm
  rw [← h, gpMap_comp]
  rfl

variable {P} in
/-- ★★★**不変量である** —— 遷移 `c` で代表元を取り替えても値が変わらない。

★★`Div c` の項が**第 1 項と第 2 項の両方に現れて打ち消す**のが要点。 -/
theorem sliceDivGpOf_precomp {A B A' A'' : C} (a : A' ⟶ A) (ha : IsIso (P.Base a))
    (has : P.degFr a = 1) (c : A'' ⟶ A') (_hc : IsIso (P.Base c))
    (hca : IsIso (P.Base (c ≫ a))) (φ : A' ⟶ B) (hcs : P.degFr c = 1) :
    sliceDivGpOf (c ≫ a) hca (c ≫ φ) = sliceDivGpOf a ha φ := by
  haveI := ha
  haveI := hca
  -- ★分子を `Φ.map (Base c)` の像としてまとめる
  have hnum : toGp _ (P.Div (c ≫ φ))
        - ((P.degFr (c ≫ φ) : ℕ+) : ℕ) • toGp _ (P.Div (c ≫ a))
      = gpMap _ (Φ.map (P.Base c))
        (toGp _ (P.Div φ) - ((P.degFr φ : ℕ+) : ℕ) • toGp _ (P.Div a)) := by
    rw [P.Div_comp, P.Div_comp, P.degFr_comp, hcs, mul_one, has, PNat.one_coe, one_smul,
      toGp_add, toGp_add, toGp_nsmul, map_sub, gpMap_toGp, map_nsmul, gpMap_toGp]
    simp only [smul_add]
    abel
  have hkey : inv (P.Base (c ≫ a)) ≫ P.Base c = inv (P.Base a) := by
    rw [IsIso.inv_comp_eq, P.Base_comp, Category.assoc, IsIso.hom_inv_id, Category.comp_id]
  rw [sliceDivGpOf_eq, sliceDivGpOf_eq, hnum, gpMap_phi_comp]
  exact congrArg (fun t => gpMap _ (Φ.map t)
    (toGp _ (P.Div φ) - ((P.degFr φ : ℕ+) : ℕ) • toGp _ (P.Div a))) hkey

/-! ## ★3. 残り —— 余錐と関手

★★不変性(上)があるので、**余錐 `biratDivGpCocone` はそのまま組める**。
★その先の `map_comp`(合成則)が `𝔽_{Φ^gp}` への関手を与え、
`Φ^birat` はその像として定まる。
★★**`map_comp` は `Prop32.lean` の `pfDiv_comp` と同じ形の計算**になる見込みで、
そこでの測定(「目標の内側を書き換えようとしない」「`congrArg₂` でまとめて外側で `rw`」)が
そのまま効くはずである。
-/

end ABC3.Found.FrdI
