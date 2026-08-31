/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.InvertibleIndex
import ABC3.Found.Arakelov.ArchDegSmul
import ABC3.Found.Arakelov.DegArith
import ABC3.Meta.Claim

/-!
# **段 A〜E がすべて閉じた** —— `deg_F` は切断に依らず、加法的である（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

## ★★★★★★★★★★これは何か

`§9-779` で **前層の大域切断 `Γ(L)` そのものが可逆加群**であることが分かったので、
`deg_F` の有限部分を**層化を経由せず**そこで測る。

    `degArithPre L̄ s = ( log #(Γ(L)/Γ(X,⊤)·s) − Σ_σ log |s|_σ ) / [F:ℚ]`

★★★これで台帳 `arakelov-degF-finite-places` の**段 A〜E がすべて閉じる**:

| 段 | 主張 | 本ファイルでの名前 |
|---|---|---|
| A | アフィンの橋 | `invertible_gammaPre`（§9-779） |
| B | 有限部分 | `degFinPre`（`finite_quotient_span_invertible`） |
| C | アルキメデス部分 | `archDeg`（在庫） |
| D | ★**切断の取り方に依らない** | `degArithPre_congr` |
| E | ★★**加法性** | `degArithPre_mul` |

## ★★★★★段 D の機構（積公式）

    有限側     `+ log |N(c)|`          （`degFinPre_smul` ＋ `card_quotient_span_eq_normAbs`）
    アルキメデス側 `− log |N(c)|/[F:ℚ]`  （`archDeg_smul`、§9-775）
    ★和        `0`

★★任意の二つの非零切断については、可逆加群の**階数 1 性**
（`exists_common_smul_invertible`）が `a·s = b·t` を与えるので、
そこに上の不変性を 2 回当てるだけである。

## ★★★★★★段 E の機構

    有限側     `card_quotient_tmul`（§9-778）——`Γ(L ⊗ M) = Γ(L) ⊗ Γ(M)` が `rfl` なので直接当たる
    アルキメデス側 `archDeg_mul`（在庫）

## ★残っている段（明示）

★★`degArithPre` はまだ**等長同型類の関数になっていない**——
`APicM` の上に降ろすには `Isometric L M → degArithPre L s = degArithPre M (φ s)` が要る。
★★★さらに `X(ℚ̄)` の上に降ろすには底変換 `deg_K(L|_K) = deg_F(L)` の計量版が要る。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite NumberField

/-! ## ★★★★`Γ(Spec R, ⊤)` と `R` の同一視 -/

/-- ★`Γ(Spec R, ⊤) ≃+* R`。 -/
noncomputable def gammaSpecRingEquiv (R : CommRingCat.{0}) :
    (Γ(Spec R, (⊤ : (Spec R).Opens)) : Type) ≃+* (R : Type) :=
  (Scheme.ΓSpecIso R).commRingCatIsoToRingEquiv

/-- ★`Γ-Spec` 同型は非零を非零に送る。 -/
theorem gammaSpec_hom_ne_zero (R : CommRingCat.{0})
    (c : (Γ(Spec R, (⊤ : (Spec R).Opens)) : Type)) (hc : c ≠ 0) :
    (Scheme.ΓSpecIso R).hom.hom c ≠ 0 := by
  intro h
  apply hc
  apply (gammaSpecRingEquiv R).injective
  show (Scheme.ΓSpecIso R).hom.hom c = (gammaSpecRingEquiv R) 0
  rw [map_zero]
  exact h

/-- ★★**主イデアルによる商の位数は `Γ-Spec` 同型で移る**。 -/
theorem card_quotient_gammaSpec (R : CommRingCat.{0})
    (c : (Γ(Spec R, (⊤ : (Spec R).Opens)) : Type)) :
    Nat.card ((Γ(Spec R, (⊤ : (Spec R).Opens)) : Type) ⧸ (Ideal.span {c}))
      = Nat.card ((R : Type)
          ⧸ (Ideal.span {(Scheme.ΓSpecIso R).hom.hom c} : Ideal (R : Type))) := by
  refine Nat.card_congr (Ideal.quotientEquiv _ _ (gammaSpecRingEquiv R) ?_).toEquiv
  rw [Ideal.map_span]
  simp
  rfl

/-- ★★`𝓞_F` の非零元による `Γ(Spec 𝓞_F, ⊤)` の商は有限である。 -/
theorem finite_quotient_gammaSpec (F : Type) [Field F] [NumberField F]
    (c : (Γ(Spec (CommRingCat.of (𝓞 F)), (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens)) : Type))
    (hc : c ≠ 0) :
    Finite ((Γ(Spec (CommRingCat.of (𝓞 F)), (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens)) : Type)
      ⧸ (Ideal.span {c})) := by
  haveI := finite_quotient_span F ((Scheme.ΓSpecIso (CommRingCat.of (𝓞 F))).hom.hom c)
    (gammaSpec_hom_ne_zero _ c hc)
  have hmap : (Ideal.span {(Scheme.ΓSpecIso (CommRingCat.of (𝓞 F))).hom.hom c}
      : Ideal ((CommRingCat.of (𝓞 F)) : Type))
      = Ideal.map (gammaSpecRingEquiv (CommRingCat.of (𝓞 F))) (Ideal.span {c}) := by
    rw [Ideal.map_span]
    simp
    rfl
  exact Finite.of_equiv _
    (Ideal.quotientEquiv _ _ (gammaSpecRingEquiv (CommRingCat.of (𝓞 F))) hmap).toEquiv.symm

/-! ## ★★★★★★★★有限部分（前層の側） -/

/-- ★★★★★★**`deg_F` の有限部分**（前層の大域切断で測る）。

    `degFinPre L s = log #(Γ(L) / Γ(X,⊤)·s)` -/
noncomputable def degFinPre {X : Scheme.{0}} (L : AInv X)
    (s : (L.carrier.sheaf.obj (op ⊤) : Type)) : ℝ :=
  Real.log (Nat.card ((L.carrier.sheaf.obj (op ⊤) : Type)
    ⧸ ((Γ(X, (⊤ : X.Opens)) : Type) ∙ s)))

/-- ★★★★★**有限部分の切断依存性**: `degFinPre (c·s) = degFinPre s + log #(Γ/(c))`。 -/
theorem degFinPre_smul {X : Scheme.{0}} [IsDomain (Γ(X, (⊤ : X.Opens)) : Type)] (L : AInv X)
    (s : (L.carrier.sheaf.obj (op ⊤) : Type)) (hs : s ≠ 0)
    (c : (Γ(X, (⊤ : X.Opens)) : Type)) (hc : c ≠ 0)
    (hfin : ∀ r : (Γ(X, (⊤ : X.Opens)) : Type), r ≠ 0 →
      Finite ((Γ(X, (⊤ : X.Opens)) : Type)
        ⧸ (Ideal.span {r} : Ideal (Γ(X, (⊤ : X.Opens)) : Type)))) :
    degFinPre L (c • s)
      = degFinPre L s
        + Real.log (Nat.card ((Γ(X, (⊤ : X.Opens)) : Type)
            ⧸ (Ideal.span {c} : Ideal (Γ(X, (⊤ : X.Opens)) : Type)))) := by
  haveI := invertible_gammaPre L
  haveI := finite_quotient_span_invertible (Γ(X, (⊤ : X.Opens)) : Type)
    (L.carrier.sheaf.obj (op ⊤) : Type) hfin s hs
  haveI := hfin c hc
  have h1 : 0 < Nat.card ((L.carrier.sheaf.obj (op ⊤) : Type)
      ⧸ ((Γ(X, (⊤ : X.Opens)) : Type) ∙ s)) := Nat.card_pos
  have h2 : 0 < Nat.card ((Γ(X, (⊤ : X.Opens)) : Type)
      ⧸ (Ideal.span {c} : Ideal (Γ(X, (⊤ : X.Opens)) : Type))) := Nat.card_pos
  show Real.log (Nat.card (_ ⧸ ((Γ(X, (⊤ : X.Opens)) : Type) ∙ (c • s)))) = _
  rw [card_quotient_smul_invertible (Γ(X, (⊤ : X.Opens)) : Type)
    (L.carrier.sheaf.obj (op ⊤) : Type) s hs c, Nat.cast_mul,
    Real.log_mul (by exact_mod_cast h1.ne') (by exact_mod_cast h2.ne')]
  rfl

/-- ★★★★★★★★★★**有限部分は加法的である**（段 E の有限側）。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★`Γ(L ⊗ M) = Γ(L) ⊗ Γ(M)` は `rfl`（前層のテンソルは各開集合ごと）なので、
§9-778 の `card_quotient_tmul` が**直接**当たる。 -/
theorem degFinPre_mul {X : Scheme.{0}} [IsDomain (Γ(X, (⊤ : X.Opens)) : Type)] (L M : AInv X)
    (s : (L.carrier.sheaf.obj (op ⊤) : Type)) (t : (M.carrier.sheaf.obj (op ⊤) : Type))
    (hs : s ≠ 0) (ht : t ≠ 0)
    (hfin : ∀ r : (Γ(X, (⊤ : X.Opens)) : Type), r ≠ 0 →
      Finite ((Γ(X, (⊤ : X.Opens)) : Type)
        ⧸ (Ideal.span {r} : Ideal (Γ(X, (⊤ : X.Opens)) : Type)))) :
    degFinPre (L.mul M) (s ⊗ₜ[(Γ(X, (⊤ : X.Opens)) : Type)] t)
      = degFinPre L s + degFinPre M t := by
  haveI := invertible_gammaPre L
  haveI := invertible_gammaPre M
  haveI hA := finite_quotient_span_invertible (Γ(X, (⊤ : X.Opens)) : Type)
    (L.carrier.sheaf.obj (op ⊤) : Type) hfin s hs
  haveI hB := finite_quotient_span_invertible (Γ(X, (⊤ : X.Opens)) : Type)
    (M.carrier.sheaf.obj (op ⊤) : Type) hfin t ht
  have h1 : 0 < Nat.card ((L.carrier.sheaf.obj (op ⊤) : Type)
      ⧸ ((Γ(X, (⊤ : X.Opens)) : Type) ∙ s)) := Nat.card_pos
  have h2 : 0 < Nat.card ((M.carrier.sheaf.obj (op ⊤) : Type)
      ⧸ ((Γ(X, (⊤ : X.Opens)) : Type) ∙ t)) := Nat.card_pos
  show Real.log (Nat.card (_ ⧸ ((Γ(X, (⊤ : X.Opens)) : Type) ∙
    (s ⊗ₜ[(Γ(X, (⊤ : X.Opens)) : Type)] t)))) = _
  rw [card_quotient_tmul (Γ(X, (⊤ : X.Opens)) : Type)
    (L.carrier.sheaf.obj (op ⊤) : Type) (M.carrier.sheaf.obj (op ⊤) : Type) s hs t,
    Nat.cast_mul, Real.log_mul (by exact_mod_cast h1.ne') (by exact_mod_cast h2.ne')]
  rfl

/-! ## ★★★★★★★★★★`deg_F`（前層の側）と段 D・E -/

/-- ★★★★★★★★★★**`deg_F(L̄)`**（前層の大域切断で測った形）。

    `degArithPre L̄ s = ( log #(Γ(L)/Γ(X,⊤)·s) − Σ_σ log |s|_σ ) / [F:ℚ]` -/
noncomputable def degArithPre (F : Type) [Field F] [NumberField F]
    (L : AInv (Spec (CommRingCat.of (𝓞 F))))
    (s : (L.carrier.sheaf.obj (op ⊤) : Type)) : ℝ :=
  degFinPre L s / (Module.finrank ℚ F : ℝ) + archDeg F L.carrier s

/-- ★★★★★★★★★**段 D の核**——切断の大域函数倍で不変（＝積公式）。 -/
theorem degArithPre_smul (F : Type) [Field F] [NumberField F]
    (L : AInv (Spec (CommRingCat.of (𝓞 F))))
    (c : (Γ(Spec (CommRingCat.of (𝓞 F)),
      (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens)) : Type))
    (s : (L.carrier.sheaf.obj (op ⊤) : Type)) (hc : c ≠ 0) (hs : s ≠ 0)
    (hnorm : ∀ σ : F →+* ℂ, L.carrier.norm s (embSpecPoint F σ) ≠ 0) :
    degArithPre F L (c • s) = degArithPre F L s := by
  have hu := gammaSpec_hom_ne_zero (CommRingCat.of (𝓞 F)) c hc
  have hfin := degFinPre_smul L s hs c hc (fun r hr => finite_quotient_gammaSpec F r hr)
  rw [card_quotient_gammaSpec (CommRingCat.of (𝓞 F)) c,
    card_quotient_span_eq_normAbs ((CommRingCat.of (𝓞 F)) : Type)
      ((Scheme.ΓSpecIso (CommRingCat.of (𝓞 F))).hom.hom c)] at hfin
  have harch := archDeg_smul F L.carrier c s hu hnorm
  show degFinPre L (c • s) / _ + archDeg F L.carrier (c • s) = _
  rw [hfin, harch]
  show _ = degFinPre L s / _ + archDeg F L.carrier s
  ring

/-- ★★★★★★★★★★**段 D が閉じた**——`deg_F` は切断の取り方に依らない。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★可逆加群の**階数 1 性**（`exists_common_smul_invertible`）が
`a·s = b·t`（`a`, `b ≠ 0`）を与えるので、
そこに `degArithPre_smul` を 2 回当てるだけである。 -/
theorem degArithPre_congr (F : Type) [Field F] [NumberField F]
    (L : AInv (Spec (CommRingCat.of (𝓞 F))))
    (s t : (L.carrier.sheaf.obj (op ⊤) : Type)) (hs : s ≠ 0) (ht : t ≠ 0)
    (hns : ∀ σ : F →+* ℂ, L.carrier.norm s (embSpecPoint F σ) ≠ 0)
    (hnt : ∀ σ : F →+* ℂ, L.carrier.norm t (embSpecPoint F σ) ≠ 0) :
    degArithPre F L s = degArithPre F L t := by
  haveI := invertible_gammaPre L
  obtain ⟨a, b, ha, hb, hab⟩ := exists_common_smul_invertible
    (Γ(Spec (CommRingCat.of (𝓞 F)), (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens)) : Type)
    (L.carrier.sheaf.obj (op ⊤) : Type) s t hs ht
  have h1 := degArithPre_smul F L a s ha hs hns
  have h2 := degArithPre_smul F L b t hb ht hnt
  rw [← h1, ← h2, hab]

/-- ★★★★★★★★★★**段 E が閉じた**——`deg_F` は加法的である。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

    `deg_F(L̄ ⊗ M̄) = deg_F(L̄) + deg_F(M̄)`

★有限側は `degFinPre_mul`（§9-778 の `card_quotient_tmul`）、
アルキメデス側は在庫の `archDeg_mul`（`|s ⊗ t| = |s|·|t|`）である。 -/
theorem degArithPre_mul (F : Type) [Field F] [NumberField F]
    (L M : AInv (Spec (CommRingCat.of (𝓞 F))))
    (s : (L.carrier.sheaf.obj (op ⊤) : Type))
    (t : (M.carrier.sheaf.obj (op ⊤) : Type)) (hs : s ≠ 0) (ht : t ≠ 0)
    (hns : ∀ σ : F →+* ℂ, L.carrier.norm s (embSpecPoint F σ) ≠ 0)
    (hnt : ∀ σ : F →+* ℂ, M.carrier.norm t (embSpecPoint F σ) ≠ 0) :
    degArithPre F (L.mul M)
        (s ⊗ₜ[(Γ(Spec (CommRingCat.of (𝓞 F)),
          (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens)) : Type)] t)
      = degArithPre F L s + degArithPre F M t := by
  have hfin := degFinPre_mul L M s t hs ht (fun r hr => finite_quotient_gammaSpec F r hr)
  have harch : archDeg F (L.mul M).carrier
      (s ⊗ₜ[(Γ(Spec (CommRingCat.of (𝓞 F)),
        (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens)) : Type)] t)
      = archDeg F L.carrier s + archDeg F M.carrier t :=
    archDeg_mul F L.carrier M.carrier s t hns hnt
  show degFinPre (L.mul M) _ / _ + archDeg F (L.mul M).carrier _ = _
  rw [hfin, harch]
  show _ = degFinPre L s / _ + archDeg F L.carrier s + (degFinPre M t / _ + archDeg F M.carrier t)
  ring

/-! ### ★出典の紐付け(`.src`) -/

def degArithPre.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(deg_F の古典的な定義——前層の大域切断で測る形)",
    sectionId := "genell-def-1-1-ii" }

def degArithPre_congr.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(段 D——deg_F は切断の取り方に依らない)",
    sectionId := "genell-def-1-1-ii" }

def degArithPre_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(段 E——deg_F は加法的である)",
    sectionId := "genell-def-1-1-ii" }

def degArithPre_mul.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "card_quotient_tmul(指数はテンソル積で掛け算になる、§9-778)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.card_quotient_tmul") 4,
    .citation "[ABC3]" "invertible_gammaPre(前層の大域切断は可逆、§9-779)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.invertible_gammaPre") 4,
    .citation "[ABC3]" "archDeg_mul(アルキメデス側の加法性、在庫)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.archDeg_mul") 4,
    .implicitStep
      ("★★残っている段: degArithPre はまだ**等長同型類の関数になっていない**。" ++
       "APicM の上に降ろすには Isometric L M → degArithPre L s = degArithPre M (φ s) が要る。" ++
       "★★★さらに X(ℚ̄) の上に降ろすには底変換 deg_K(L|_K) = deg_F(L) の計量版が要る") 4 ]

end ABC3.Found.Arakelov
