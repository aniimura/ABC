/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.SemistableFin
import ABC3.Found.GaloisRep.HtFaltBounds
import ABC3.Found.GaloisRep.VeluNormalized
import ABC3.Found.GaloisRep.DegInfBaseChange
import ABC3.Meta.Claim
import ABC3.Found.GaloisRep.Lemma35Concrete.Lemma32

/-!
# Lemma35Concrete —— `[GenEll] Lemma 3.5` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain NumberField WeierstrassCurve ABC3.Found.GenEll
open scoped Classical
variable {L : Type} [Field L] [NumberField L]

/-- ★★★★★**悪い還元の素点だけで確かめれば足りる形**。

★良い還元の素点では `v_p(Δ_min) = 0`（両方）なので関係は自動的に成り立つ。
★★`Lemma 3.5` の設定では `E` と `E′` は同種なので悪い素点の集合 `S` は共通である。 -/
theorem degInfOf_eq_of_local' (E E' : WeierstrassCurve L) (l : ℕ)
    (S : Set (HeightOneSpectrum (𝓞 L)))
    (hbad : ∀ p ∈ S, minDeltaExp p E' = l * minDeltaExp p E)
    (hgoodE : ∀ p ∉ S, minDeltaExp p E = 0)
    (hgoodE' : ∀ p ∉ S, minDeltaExp p E' = 0) :
    degInfOf L E' = (l : ℝ) * degInfOf L E := by
  refine degInfOf_eq_of_local E E' l (fun p => ?_)
  by_cases hp : p ∈ S
  · exact hbad p hp
  · rw [hgoodE p hp, hgoodE' p hp, mul_zero]

/-! ## ★★★★★★★★★★★★★★`Lemma 3.5` —— 穴 1 つを残して -/

/-- ★★★★★★★★★★★★★★**[GenEll] Lemma 3.5**（Global Rank One Subgroups of
`l`-Torsion）——★**同種写像の高さ評価 `hfalt` だけを仮説として受けた形**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

    `(1/(12(1+ϵ)))·l·deg_∞(E) ≤ ht^Falt(E) + 2·log(l) + C`

★★`C` は `ϵ` にしか依らない——`L`・`E`・`E′`・`l` に依らない（原文が明記している）。

★入力 3 つのうち 2 つは**証明済み**である:

| 入力 | 状態 |
|---|---|
| `hdeg`（`deg∞(E′) = l·deg∞(E)`） | ★`Lemma 3.2, (ii)`（§9-1011）＋本ファイルの `degInfOf_eq_of_local` |
| `Proposition 3.4` を `E′` に当てる | ★★**済**（§9-1009、証明の中で使っている） |
| `hfalt`（`ht^Falt(E′) ≤ ht^Falt(E) + 2log l`） | ☆**残る唯一の入力** |

☆`hfalt` は本プロジェクトのどこにも証明が無い外部引用（[FC] Ch. I, Prop 2.7 ＋
アルキメデスの (1,1)-形式）なので、**項目全体の `.src` は置かない**。 -/
theorem lemma_3_5_of_isogeny_estimate (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L] (E E' : WeierstrassCurve L)
      [E.IsElliptic] [E'.IsElliptic] (l : ℕ),
      (∀ p, SemistableAt p E') →
      degInfOf L E' = (l : ℝ) * degInfOf L E →
      htFaltOf L E' ≤ htFaltOf L E + 2 * Real.log l →
      (1 / (12 * (1 + eps))) * (l : ℝ) * degInfOf L E
        ≤ htFaltOf L E + 2 * Real.log l + C := by
  obtain ⟨C, hC⟩ := prop_3_4_chain_semistable eps heps
  have hpos : (0:ℝ) < 12 * (1 + eps) := by linarith
  refine ⟨C / (12 * (1 + eps)), fun L _ _ E E' _ _ l hss hdeg hfalt => ?_⟩
  obtain ⟨h1, h2, _⟩ := hC L E' hss
  have hkey : degInfOf L E' ≤ 12 * (1 + eps) * htFaltOf L E' + C := le_trans h1 h2
  rw [hdeg] at hkey
  have hstep : (l : ℝ) * degInfOf L E
      ≤ 12 * (1 + eps) * (htFaltOf L E + 2 * Real.log l) + C := by
    refine le_trans hkey ?_
    have := mul_le_mul_of_nonneg_left hfalt hpos.le
    linarith
  rw [show (1 / (12 * (1 + eps))) * (l : ℝ) * degInfOf L E
      = ((l : ℝ) * degInfOf L E) / (12 * (1 + eps)) by ring,
    div_le_iff₀ hpos]
  have hC' : C / (12 * (1 + eps)) * (12 * (1 + eps)) = C := by field_simp
  nlinarith [hstep, hC']

/-! ## ★★★★★★★★★★★★★★★★★★残っていた入力 `hfalt` を埋める -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] Lemma 3.5 —— `hfalt` を外した形**

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

    `(1/(12(1+ϵ)))·l·deg_∞(E) ≤ ht^Falt(E) + 2·log(l) + C`

★★★★★★★☆**かつて「残る唯一の入力」だった `hfalt`
（`ht^Falt(E′) ≤ ht^Falt(E) + 2·log(l)`、[FC] Ch. I, Prop 2.7 ＋
アルキメデスの (1,1)-形式）が、`htFalt_veluQuotientFull_le`（`§9-1146`、第 704）で
埋まった。** これで本主張の入力から**未証明の外部引用が消えた**。

★受けている仮説はすべて**幾何のデータ**である:

| 仮説 | 意味 | 出どころ |
|---|---|---|
| `hQ` | `Q` は位数ちょうど `l` | 原文の `l`-捻れ部分群 |
| `hE'` | `E′ = E/⟨Q⟩`（Vélu） | 定義 |
| `P`・`hΔ`・`hPC` | 各 `σ` での一意化 | `exists_periodPair_of_isElliptic`（第 348、無条件） |
| `hmin`・`hint` | `E` が極小・`E′` が整 | モデルの取り方 |
| `hss` | `E′` は半安定 | 原文の設定 |
| `hdeg` | `deg∞(E′) = l·deg∞(E)` | `Lemma 3.2, (ii)` ＋ `degInfOf_eq_of_local` |

☆`hdeg` はまだ仮説である（`Lemma 3.2, (ii)` の大域化に還元済み）ので、
**項目全体の `.src` はまだ置かない**（`check.mjs` B6）。 -/
theorem lemma_3_5_velu (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L] (E E' : WeierstrassCurve L)
      [E.IsElliptic] [E'.IsElliptic] (l : ℕ), 0 < l →
      ∀ Q : E.toAffine.Point, addOrderOf Q = l →
      E' = veluQuotientFull E (((Finset.range l).erase 0).image
          (fun k : ℕ => pointCoords (k • Q))) →
      ∀ (P : (L →+* ℂ) → PeriodPair) (Cv : (L →+* ℂ) → VariableChange ℂ),
      (∀ σ, latticeDisc (P σ) ≠ 0) →
      (∀ σ, Cv σ • (E.map σ) = latticeCurve (P σ)) →
      (∀ σ : L →+* ℂ, (E.map σ).IsElliptic) →
      (∀ σ : L →+* ℂ, (Cv σ • (E.map σ)).IsElliptic) →
      (∀ p : HeightOneSpectrum (𝓞 L), neronExp p E = 0) →
      (∀ p : HeightOneSpectrum (𝓞 L), E'.IsIntegral (primeSubring p)) →
      (∀ p, SemistableAt p E') →
      degInfOf L E' = (l : ℝ) * degInfOf L E →
      (1 / (12 * (1 + eps))) * (l : ℝ) * degInfOf L E
        ≤ htFaltOf L E + 2 * Real.log l + C := by
  obtain ⟨C, hC⟩ := lemma_3_5_of_isogeny_estimate eps heps
  refine ⟨C, fun L _ _ E E' _ _ l hl Q hQ hE' P Cv hΔ hPC hell1 hell2 hmin hint
    hss hdeg => ?_⟩
  exact hC L E E' l hss hdeg
    (htFalt_veluQuotientFull_le E E' l hl Q hQ hE' P Cv hΔ hPC hell1 hell2 hmin hint)

def lemma_3_5_velu.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(hfalt を外した形——未証明の外部引用が消えた。残るは hdeg のみ)",
    sectionId := "genell-lemma-3-5" }

def lemma_3_5_velu.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "htFalt_veluQuotientFull_le(hfalt を埋めた段、§9-1146)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.htFalt_veluQuotientFull_le") 4,
    .citation "[ABC3]" "lemma_3_5_of_isogeny_estimate(組み立て)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.lemma_3_5_of_isogeny_estimate") 3,
    .implicitStep
      ("☆残るのは hdeg(deg∞(E′) = l·deg∞(E))を Lemma 3.2, (ii) の局所版から " ++
       "degInfOf_eq_of_local で大域化して外すことである") 3 ]

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] Lemma 3.5 —— 残る入力は局所の `v_p(Δ_min(E′)) = l·v_p(Δ_min(E))` だけ**

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★★★★★☆第 705 で `hfalt` が消えたので、**残っているのは
`hloc`（各素点での `Δ_min` の関係）だけ**である。これは `Lemma 3.2, (ii)`
（`Found/GenEll/Lemma32.lean`、`vAdd v (q^l) = l·vAdd v q`）が与えるもので、
Tate 曲線の母数 `q_E` と `Δ_min` の対応を挟めば出る。

★★☆**外部文献への未証明の依存はもう無い**——残るのは本プロジェクト内の接続だけである。 -/
theorem lemma_3_5_velu_local (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L] (E E' : WeierstrassCurve L)
      [E.IsElliptic] [E'.IsElliptic] (l : ℕ), 0 < l →
      ∀ Q : E.toAffine.Point, addOrderOf Q = l →
      E' = veluQuotientFull E (((Finset.range l).erase 0).image
          (fun k : ℕ => pointCoords (k • Q))) →
      ∀ (P : (L →+* ℂ) → PeriodPair) (Cv : (L →+* ℂ) → VariableChange ℂ),
      (∀ σ, latticeDisc (P σ) ≠ 0) →
      (∀ σ, Cv σ • (E.map σ) = latticeCurve (P σ)) →
      (∀ σ : L →+* ℂ, (E.map σ).IsElliptic) →
      (∀ σ : L →+* ℂ, (Cv σ • (E.map σ)).IsElliptic) →
      (∀ p : HeightOneSpectrum (𝓞 L), neronExp p E = 0) →
      (∀ p : HeightOneSpectrum (𝓞 L), E'.IsIntegral (primeSubring p)) →
      (∀ p, SemistableAt p E') →
      (∀ p : HeightOneSpectrum (𝓞 L), minDeltaExp p E' = l * minDeltaExp p E) →
      (1 / (12 * (1 + eps))) * (l : ℝ) * degInfOf L E
        ≤ htFaltOf L E + 2 * Real.log l + C := by
  obtain ⟨C, hC⟩ := lemma_3_5_velu eps heps
  refine ⟨C, fun L _ _ E E' _ _ l hl Q hQ hE' P Cv hΔ hPC hell1 hell2 hmin hint
    hss hloc => ?_⟩
  exact hC L E E' l hl Q hQ hE' P Cv hΔ hPC hell1 hell2 hmin hint hss
    (degInfOf_eq_of_local E E' l hloc)

def lemma_3_5_velu_local.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(残る入力は局所の v_p(Δ_min(E′)) = l·v_p(Δ_min(E)) だけ)",
    sectionId := "genell-lemma-3-5" }

def lemma_3_5_velu_local.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "lemma_3_5_velu(hfalt を外した形、§9-1147)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.lemma_3_5_velu") 2,
    .citation "[ABC3]" "degInfOf_eq_of_local(局所から大域へ)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.degInfOf_eq_of_local") 2,
    .otherPaper "[GenEll]" "Lemma 3.2, (ii)(vAdd v (q^l) = l·vAdd v q——局所の Δ_min の関係)" 15 ]

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] Lemma 3.5 —— 残る入力は**悪い還元の素点だけ**の `Δ_min` の関係**

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★良い還元の素点では `v_p(Δ_min) = 0`（両方）なので関係は自動である。
★★★★★★★★★☆したがって残っているのは

    `∀ p ∈ S, v_p(Δ_min(E′)) = l·v_p(Δ_min(E))`（`S` は悪い還元の素点）

だけ——これは `Lemma 3.2, (ii)`（Tate 曲線の母数 `q_{E′} = q_E^l`）そのものである。 -/
theorem lemma_3_5_velu_bad (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L] (E E' : WeierstrassCurve L)
      [E.IsElliptic] [E'.IsElliptic] (l : ℕ), 0 < l →
      ∀ Q : E.toAffine.Point, addOrderOf Q = l →
      E' = veluQuotientFull E (((Finset.range l).erase 0).image
          (fun k : ℕ => pointCoords (k • Q))) →
      ∀ (P : (L →+* ℂ) → PeriodPair) (Cv : (L →+* ℂ) → VariableChange ℂ),
      (∀ σ, latticeDisc (P σ) ≠ 0) →
      (∀ σ, Cv σ • (E.map σ) = latticeCurve (P σ)) →
      (∀ σ : L →+* ℂ, (E.map σ).IsElliptic) →
      (∀ σ : L →+* ℂ, (Cv σ • (E.map σ)).IsElliptic) →
      (∀ p : HeightOneSpectrum (𝓞 L), neronExp p E = 0) →
      (∀ p : HeightOneSpectrum (𝓞 L), E'.IsIntegral (primeSubring p)) →
      (∀ p, SemistableAt p E') →
      ∀ S : Set (HeightOneSpectrum (𝓞 L)),
      (∀ p ∈ S, minDeltaExp p E' = l * minDeltaExp p E) →
      (∀ p ∉ S, minDeltaExp p E = 0) →
      (∀ p ∉ S, minDeltaExp p E' = 0) →
      (1 / (12 * (1 + eps))) * (l : ℝ) * degInfOf L E
        ≤ htFaltOf L E + 2 * Real.log l + C := by
  obtain ⟨C, hC⟩ := lemma_3_5_velu eps heps
  refine ⟨C, fun L _ _ E E' _ _ l hl Q hQ hE' P Cv hΔ hPC hell1 hell2 hmin hint
    hss S hbad hgoodE hgoodE' => ?_⟩
  exact hC L E E' l hl Q hQ hE' P Cv hΔ hPC hell1 hell2 hmin hint hss
    (degInfOf_eq_of_local' E E' l S hbad hgoodE hgoodE')

def lemma_3_5_velu_bad.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(残る入力は悪い還元の素点だけの Δ_min の関係)",
    sectionId := "genell-lemma-3-5" }

def lemma_3_5_velu_bad.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "lemma_3_5_velu(hfalt を外した形、§9-1147)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.lemma_3_5_velu") 2,
    .citation "[ABC3]" "degInfOf_eq_of_local'(悪い素点だけで確かめれば足りる形)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.degInfOf_eq_of_local'") 1,
    .otherPaper "[GenEll]" "Lemma 3.2, (ii)(q_{E′} = q_E^l——悪い素点での Δ_min の関係)" 15 ]

/-! ## ★★★★★★★★★★★★★★★★★★残る入力は **`j` の付値だけ**である -/

/-- ★★★★★★★★★★★★**半安定なら `v_p(Δ_min)` の関係は `v_p(j)` の関係から出る**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★半安定なら `v_p(Δ_min) = max(0, −v_p(j))`（`§9-1165`、第 738）なので、
`v_p(j′) = l·v_p(j)` から直ちに `v_p(Δ_min(E′)) = l·v_p(Δ_min(E))` が出る。 -/
theorem minDeltaExp_of_jExp_mul {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    (hss : SemistableAt p E) (hss' : SemistableAt p E') (l : ℕ)
    (hbad : jExp p E < 0 → jExp p E' = (l : ℤ) * jExp p E)
    (hgood : 0 ≤ jExp p E → 0 ≤ jExp p E') :
    minDeltaExp p E' = (l : ℤ) * minDeltaExp p E := by
  rw [minDeltaExp_eq_maxJ_of_semistable p E' hss', minDeltaExp_eq_maxJ_of_semistable p E hss]
  rcases lt_or_ge (jExp p E) 0 with h | h
  · rw [hbad h, ← mul_neg, max_zero_mul _ (by positivity)]
  · have hg := hgood h
    have h1 : max 0 (-jExp p E') = 0 := max_eq_left (by omega)
    have h2 : max 0 (-jExp p E) = 0 := max_eq_left (by omega)
    rw [h1, h2, mul_zero]

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`Lemma 3.5`——残る入力は**悪い素点での** `v_p(j′) = l·v_p(j)` **だけ**である**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★これが本項目の**最終形**である。`hfalt`（Faltings 高さの評価）は
`§9-1164`、第 704 で内製し、`hdeg`（`deg∞(E′) = l·deg∞(E)`）は
半安定性と合わせて**`j` の付値の関係 1 本**に落ちた。

☆残る 1 本は `Skeleton/GenEll/TateIsogeny.lean` の `tateModel_of_quot_mu`
（`E_q/μ_l` は母数 `q^l` の Tate 曲線）である
——Tate 曲線では `v_p(j) = −v_p(q)` なので、それがそのまま `v_p(j′) = l·v_p(j)` を与える。

☆良い還元の素点（`v_p(j) ≥ 0`）では `v_p(j′) ≥ 0` だけでよい
——半安定なら `v_p(Δ_min) = max(0, −v_p(j))` なので両辺とも `0` になる。
★ここを `v_p(j′) = l·v_p(j)` と書くのは**強すぎ**である（`v_p(j) > 0` のとき假になりうる）。 -/
theorem lemma_3_5_velu_j (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L] (E E' : WeierstrassCurve L)
      [E.IsElliptic] [E'.IsElliptic] (l : ℕ), 0 < l →
      ∀ Q : E.toAffine.Point, addOrderOf Q = l →
      E' = veluQuotientFull E (((Finset.range l).erase 0).image
          (fun k : ℕ => pointCoords (k • Q))) →
      ∀ (P : (L →+* ℂ) → PeriodPair) (Cv : (L →+* ℂ) → VariableChange ℂ),
      (∀ σ, latticeDisc (P σ) ≠ 0) →
      (∀ σ, Cv σ • (E.map σ) = latticeCurve (P σ)) →
      (∀ σ : L →+* ℂ, (E.map σ).IsElliptic) →
      (∀ σ : L →+* ℂ, (Cv σ • (E.map σ)).IsElliptic) →
      (∀ p : HeightOneSpectrum (𝓞 L), neronExp p E = 0) →
      (∀ p : HeightOneSpectrum (𝓞 L), E'.IsIntegral (primeSubring p)) →
      (∀ p, SemistableAt p E) →
      (∀ p, SemistableAt p E') →
      (∀ p : HeightOneSpectrum (𝓞 L), jExp p E < 0 → jExp p E' = (l : ℤ) * jExp p E) →
      (∀ p : HeightOneSpectrum (𝓞 L), 0 ≤ jExp p E → 0 ≤ jExp p E') →
      (1 / (12 * (1 + eps))) * (l : ℝ) * degInfOf L E
        ≤ htFaltOf L E + 2 * Real.log l + C := by
  obtain ⟨C, hC⟩ := lemma_3_5_velu eps heps
  refine ⟨C, fun L _ _ E E' _ _ l hl Q hQ hE' P Cv hΔ hPC hell1 hell2 hmin hint
    hssE hssE' hbad hgood => ?_⟩
  refine hC L E E' l hl Q hQ hE' P Cv hΔ hPC hell1 hell2 hmin hint hssE'
    (degInfOf_eq_of_local E E' l (fun p => ?_))
  exact minDeltaExp_of_jExp_mul p E E' (hssE p) (hssE' p) l (hbad p) (hgood p)

def minDeltaExp_of_jExp_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(半安定なら Δ_min の関係は j の付値の関係から出る。★無条件)",
    sectionId := "genell-lemma-3-5" }

def lemma_3_5_velu_j.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(残る入力は悪い素点での v_p(j′) = l·v_p(j) だけ——最終形)",
    sectionId := "genell-lemma-3-5" }

def lemma_3_5_velu_j.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "tateModel_of_quot_mu(E_q/μ_l は母数 q^l の Tate 曲線)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.tateModel_of_quot_mu") 15,
    .implicitStep
      ("☆Tate 曲線では v_p(j) = −v_p(q) なので、母数が q^l であることが" ++
       "そのまま v_p(j′) = l·v_p(j) を与える") 3 ]

def lemma_3_5_of_isogeny_estimate.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(曲線の水準で——同種写像の高さ評価だけを仮説として受けた形)",
    sectionId := "genell-lemma-3-5" }

def lemma_3_5_of_isogeny_estimate.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "prop_3_4_chain_semistable(Prop 3.4 を E′ に当てる段、§9-1009)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.prop_3_4_chain_semistable") 4,
    .citation "[ABC3]" "lemma_3_2(deg∞(E′) = l·deg∞(E) の局所版、§9-1011)"
      (.inProject "ABC3" "ABC3.Found.GenEll.lemma_3_2") 4,
    .otherPaper "[FC]"
      ("Chapter I, Proposition 2.7 ＋ 原文「since integrating a (1,1)-form over Ev " ++
       "differs from integrating over (EH)v by a factor of l」" ++
       "——**ht^Falt(E′) ≤ ht^Falt(E) + 2·log(l)**。" ++
       "☆★★★★★★★★★★**Lemma 3.5 に残る唯一の入力**であり、" ++
       "本ファイルは仮説 hfalt として**型で固定した**。" ++
       "★受けたままでは項目に数えない(check.mjs B6)。" ++
       "★★mathlib には楕円曲線の同種写像も Néron モデルも無い" ++
       "(2026-08-29 に #check で確認: WeierstrassCurve.Isogeny・" ++
       "AlgebraicGeometry.NeronModel・EllipticCurve.Point いずれも不在)") 9,
    .implicitStep
      ("★★★★★★到達点(2026-08-29、第 565): Lemma 3.5 の 3 入力のうち**2 つが証明済み**になり、" ++
       "残る 1 つを**曲線の水準の型のついた仮説**として Found/ に固定した。" ++
       "★Proposition 1.7 が Skeleton で hlow／hup を受けてから " ++
       "§9-954〜§9-974 の 21 ブロックで作ったのと同じ順序である") 8,
    .implicitStep
      ("☆`degInfOf_eq_of_local'` が要求する「E と E′ の悪い素点の集合が共通」は、" ++
       "同種な曲線が同じ還元型を持つこと(Néron–Ogg–Shafarevich)による。" ++
       "★半安定の設定では Lemma 3.2, (ii) が乗法還元の素点で関係を与え、" ++
       "良還元の素点では両方 0 である") 6 ]

end ABC3.Found.GaloisRep
