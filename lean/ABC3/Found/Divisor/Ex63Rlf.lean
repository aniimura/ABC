/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.Ex63Types

/-!
# [FrdI] Theorem 6.4, (i) —— 算術の `𝒞^rlf` は**もう 1 つの `ArithDatum`** で建つ

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.114。

原文 (FrdI p.114):
> Then the Frobenioids C, Cpf, Crlf, Cun-tr, (Cpf)un-tr are of isotropic

## ★★★★★★★これが `𝒞^rlf` を止めていた壁の**迂回路**である

`𝒞^rlf` は一般には `Φ^rlf = ℝ≥0 ⊗_ℕ Φ` （`ScT`）で作るが、そこには

* `hcharInj` —— `f : M ↪ N ⟹ id ⊗ f` が単射（節点 `rlf-flat`）
* `hdivRlf` —— `ScT ℝ≥0 Φ(A)` が divisorial（簡約性・飽和性・characteristic 型）

の 2 つが要り、★★**どちらも「`ℝ≥0` が ℕ 上平坦」という凸幾何**に帰着する
（`hcharInj` は算術ではトレース写像で迂回できたが、`hdivRlf` は迂回できない）。

### ★★★実現化した `Φ^rlf(L)` は**具体的に書ける**

`Φ(L) = Γ_L ∩ ℝ≥0[V(L)]`（`Γ_L = arithDivGroup L`、有限素点で `ℤ·log N𝔭`）を
実現化したものは、**全素点で実係数の有効因子**

    Φ^rlf(L) = ℝ≥0[V(L)] = effR (⊤ : AddSubgroup (ArithPlace L →₀ ℝ))

に他ならない。★★これは **`grp := ⊤` とした `ArithDatum`** そのものである。

したがって `𝒞^rlf` は `Example 6.3` の機械を**もう一度**、`Γ_L := ⊤` で回すだけで建つ。
`ScT` も `rlfCone` も、平坦性も要らない。

| 条件 | `arithDatumGalois` | ★`arithDatumRlf`（`grp = ⊤`） |
|---|---|---|
| `gen`（各素点に正の元） | `isGenSubgroupR_arithDivGroup` | ★自明 |
| `locMono`（局所群は離散か全体） | 有限素点は離散・無限素点は全体 | ★**全部が全体** |
| `coord`（座標ごとに閉じる） | `isCoordwiseR_arithDivGroup` | ★自明 |
| divisorial | `isDivisorial_effR` | ★同じ（**無料**） |
| perf-factorial | `isPerfFactorial_effR` | ★同じ（**無料**） |

## ★★★★逸脱の記録（CLAUDE.md の「逸脱」）

★原典は `𝒞^rlf` を「`𝒞` の実現化」として定義する。本ファイルは算術の場合に
**実現化した対象を直接 `ArithDatum` として書き下す**読み替えを行う。

★★逸脱の内容: **`Crlf(𝒞_{F̄/F}) ≅ 𝒞(Δ^rlf)` を証明していない**
（それには `ℝ≥0 ⊗_ℕ Φ(L) ≅ effR ⊤`、すなわち平坦性が要る）。
★理由: 後続（`Theorem 6.4` の (ii)(iii)(iv) と `Lemma 6.5`）が使うのは
「`Φ^rlf(L)` が全素点上の実係数である」ことだけで、
テンソル積としての表示は使わないためである。
★★`δ_A : Pic_Φ(A) ≅ ℝ`（`arithPicIso`、`ArithPicR.lean`）も
**`(ArithPlace L →₀ ℝ) ⧸ principalSpan L`** の形で書かれており、
まさにこの `Δ^rlf` の水準の主張である。
-/

namespace ABC3.Found.Divisor

open ABC3.Found.FrdI CategoryTheory

/-! ## ★1. `grp = ⊤` の `ArithDatum` -/

variable (F Kbar : Type) [Field F] [NumberField F] [Field Kbar] [Algebra F Kbar]

/-- ★★★★★★**`𝒞^rlf` の算術のデータ** —— `Φ^rlf(L) = ℝ≥0[V(L)]`
（全素点で実係数の有効因子）。

★`arithDatumGalois` との違いは `grp` を `arithDivGroup L` から `⊤` に替えた 1 点だけで、
引き戻し（`pullOf`）はそのまま使える。 -/
noncomputable def arithDatumRlf : ArithDatum.{0, 0, 0} (FinSub F Kbar)ᵒᵖ where
  primes A := ArithPlace A.unop.toIF
  pull {_ _} α := pullOf α.unop
  pull_id A x := pullOf_id A.unop x
  pull_comp {_ _ _} α β x := pullOf_comp α.unop β.unop x
  grp _ := ⊤
  pull_mem {_ _} _ {_} _ := AddSubgroup.mem_top _
  pull_nonneg {_ _} α {_} hx := by
    letI := algOfHom α.unop
    exact arithExtend_nonneg hx
  pull_inj {_ _} α := by
    letI := algOfHom α.unop
    exact arithExtend_injective
  gen _ := fun _ => ⟨1, one_pos, AddSubgroup.mem_top _⟩
  locMono _ := fun _ => Or.inr (fun _ => AddSubgroup.mem_top _)
  coord _ := fun _ _ _ => AddSubgroup.mem_top _

variable {F Kbar}

/-- ★**`Φ^rlf(L)` は全素点で実係数の有効因子**（定義的に等しい）。 -/
theorem arithDatumRlf_phi_val [IsGalois F Kbar] (A : (FinSub F Kbar)ᵒᵖ) :
    ((arithDatumRlf F Kbar).phi finSubOp_isOfFSMType).val A
      = effR (⊤ : AddSubgroup (ArithPlace A.unop.toIF →₀ ℝ)) := rfl

variable (F Kbar) in
/-- ★★★**`Φ^rlf` は divisorial**（★`ScT` の簡約性・飽和性は要らない）。 -/
theorem arithDatumRlf_isDivisorialOn [IsGalois F Kbar] :
    ((arithDatumRlf F Kbar).phi finSubOp_isOfFSMType).IsDivisorialOn :=
  (arithDatumRlf F Kbar).phi_isDivisorialOn finSubOp_isOfFSMType

variable (F Kbar) in
/-- ★★★**`Φ^rlf(L)` は perf-factorial**。 -/
theorem arithDatumRlf_isPerfFactorial [IsGalois F Kbar] (A : (FinSub F Kbar)ᵒᵖ) :
    IsPerfFactorial (((arithDatumRlf F Kbar).phi finSubOp_isOfFSMType).val A) :=
  (arithDatumRlf F Kbar).phi_isPerfFactorial finSubOp_isOfFSMType A

/-! ## ★2. 有理関数の除数 `div : B(L) → Φ^rlf(L)^gp` -/

/-- ★`Lˣ → ⊤` —— `arithDiv` をそのまま（`⊤` なので所属は自明）。 -/
noncomputable def arithDivGroupHomTop (L : Type) [Field L] [NumberField L] :
    Additive Lˣ →+ (⊤ : AddSubgroup (ArithPlace L →₀ ℝ)) :=
  AddMonoidHom.codRestrict
    (AddMonoidHom.comp ((arithDivGroup L).subtype) (arithDivGroupHom L)) _
    (fun _ => AddSubgroup.mem_top _)

variable [IsGalois F Kbar]

/-- ★★**`div : B(L) → Φ^rlf(L)^gp`**。 -/
noncomputable def divBRlf (A : (FinSub F Kbar)ᵒᵖ) :
    (bmonGalois F Kbar).val A →+
      Gp (((arithDatumRlf F Kbar).phi finSubOp_isOfFSMType).val A) :=
  ((effRGpEquiv (⊤ : AddSubgroup (ArithPlace A.unop.toIF →₀ ℝ))
      ((arithDatumRlf F Kbar).gen A)).symm.toAddMonoidHom).comp
    (arithDivGroupHomTop A.unop.toIF)

theorem phiGpHom_divBRlf (A : (FinSub F Kbar)ᵒᵖ) (x : (bmonGalois F Kbar).val A) :
    Subtype.val (phiGpHom (arithDatumRlf F Kbar) finSubOp_isOfFSMType (divBRlf A x))
      = arithDiv (Additive.toMul x) := by
  show Subtype.val ((effRGpEquiv (⊤ : AddSubgroup (ArithPlace A.unop.toIF →₀ ℝ))
      ((arithDatumRlf F Kbar).gen A))
      ((effRGpEquiv (⊤ : AddSubgroup (ArithPlace A.unop.toIF →₀ ℝ))
        ((arithDatumRlf F Kbar).gen A)).symm
        (arithDivGroupHomTop A.unop.toIF x))) = _
  rw [AddEquiv.apply_symm_apply]
  rfl

/-- ★★★**`div` は関手的** —— 中身は `arithExtend_arithDiv` 1 本（`Φ` と同じ）。 -/
theorem divBRlf_nat {A B : (FinSub F Kbar)ᵒᵖ} (f : A ⟶ B) (x : (bmonGalois F Kbar).val B) :
    divBRlf A ((bmonGalois F Kbar).map f x)
      = ((arithDatumRlf F Kbar).phi finSubOp_isOfFSMType).gpMapOn f (divBRlf B x) := by
  refine phiGpHom_injective (arithDatumRlf F Kbar) finSubOp_isOfFSMType (Subtype.ext ?_)
  rw [phiGpHom_divBRlf, phiGpHom_gpMapOn, phiGpHom_divBRlf]
  letI := algOfHom f.unop
  show arithDiv (Units.map
      (algebraMap (Opposite.unop B).toIF (Opposite.unop A).toIF).toMonoidHom (Additive.toMul x))
    = arithExtend (L := (Opposite.unop B).toIF) (arithDiv (Additive.toMul x))
  exact arithExtend_arithDiv (Additive.toMul x)

/-! ## ★3. model Frobenioid -/

variable (F Kbar)

/-- ★★★★★**`𝒞^rlf` の `ModelData`** —— `Theorem 5.2` の入力そのもの。 -/
noncomputable def rlfModelData : ModelData.{0, 0, 0} (FinSub F Kbar)ᵒᵖ :=
  (arithDatumRlf F Kbar).modelData finSubOp_isOfFSMType (bmonGalois F Kbar)
    (fun A => divBRlf A) (fun f x => divBRlf_nat f x)

/-- ★★★★★★**算術の `𝒞^rlf` は Frobenioid である**（★平坦性は要らない）。 -/
theorem rlfHyp : (rlfModelData F Kbar).Hyp :=
  (arithDatumRlf F Kbar).modelHyp finSubOp_isOfFSMType (bmonGalois F Kbar)
    (fun A => divBRlf A) (fun f x => divBRlf_nat f x)
    (bmonGalois_isGroupLike F Kbar) finSubOp_totallyEpimorphic finSubOp_isConnected

/-! ### ★出典の紐付け -/

def arithDatumRlf.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — 𝒞^rlf の算術のデータ(Φ^rlf(L) = ℝ≥0[V(L)])",
    sectionId := "frdi-thm-6-4" }

def arithDatumRlf.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("★逸脱の記録: 原典は 𝒞^rlf を「𝒞 の実現化」として定義するが、算術の場合に " ++
       "実現化した対象を直接 ArithDatum(grp = ⊤)として書き下す読み替えを行った。" ++
       "Crlf(𝒞) ≅ 𝒞(Δ^rlf) は証明していない(それには ℝ≥0 の ℕ 上の平坦性が要る)。" ++
       "後続が使うのは「Φ^rlf(L) が全素点上の実係数である」ことだけである") 114 ]

def rlfHyp.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — 算術の 𝒞^rlf は Frobenioid である(平坦性不要)",
    sectionId := "frdi-thm-6-4" }

def arithDatumRlf_isDivisorialOn.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — Φ^rlf は divisorial(ScT の簡約性は要らない)",
    sectionId := "frdi-thm-6-4" }

def arithDatumRlf_isPerfFactorial.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — Φ^rlf(L) は perf-factorial",
    sectionId := "frdi-thm-6-4" }

end ABC3.Found.Divisor
