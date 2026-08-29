/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.SemistableFin
import ABC3.Found.GaloisRep.NorthcottHtJ
import ABC3.Found.GenEll.BDClass
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★`Proposition 3.4`（`Found`、項目まるごと）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★これは何か

`Proposition 3.4` を、**半安定楕円曲線の族にわたる 1 本の主張**として取る:

    `deg_∞ ≳ ht_∞`,  `ht_∞ ≳ 12(1+ϵ)·ht^Falt`,  `12(1+ϵ)·ht^Falt ≳ (1+ϵ)·ht_∞`
    ＋ **`{ht^Falt ≤ C}` の `j` の像は有限**

★★★★★★定数はすべて **`ϵ` にしか依らない**——族にも体にも曲線にも依らない。

## ★★★量の対応（すべて `Found/` の具体的な量）

| 原文 | 本ファイル | どこ |
|---|---|---|
| `deg_∞` | `degInfOf L E = (1/d)Σ_p v_p(Δ_min)·log N(p)` | `Found/GaloisRep/DegInf.lean`（§9-319） |
| `ht^Falt` | `htFaltOf L E`（`12ht^Falt = deg∞ − archSum/d`） | `Found/GaloisRep/FaltingsWitness.lean`（§9-670） |
| `ht_∞` | `htJ L E`（`j` の素朴 Weil 高さ） | `Found/GaloisRep/HtJBound.lean`（§9-1002）★逸脱、下記 |
| `M_ell(ℚ̄)^{≤d}` の点 | `j` 不変量（`[L:ℚ] ≤ d`） | ★逸脱、下記 |

## ★★2026-08-29 の前進（明示）

| | `Skeleton/GenEll/Section3.lean` の `prop_3_4` | ★本ファイル |
|---|---|---|
| 3 本の `≲` | `EllModuliData` の欄（**posit**）から導く | ★**具体的な量で証明する** |
| 有限性 | `D.northcott`（**posit**） | ★**証明する**（Northcott） |
| `ht_∞` | 名前だけ | ★**`j` の Weil 高さとして構成**（mathlib の `Height.logHeight₁`） |
| `ht^Falt` の中身 | 名前だけ | ★**Faltings 高さの式**（`§9-670`） |

★★`Skeleton` の `prop_3_4` は界面の 4 つの公理から**組み立てるだけ**だった。
★★★本ファイルはその 4 つの公理の**中身を作る**
——`§9-1000`〜`§9-1008` の 9 ブロックがそれである。

## ★★★★★★★★★★逸脱（明示）

| 項 | 原典 | 形式化 | 理由 |
|---|---|---|---|
| `ht_∞` | `M̄_ell` の**無限遠因子に付随する高さ** | ★**`j` の素朴 Weil 高さ `h(j)`** | `M̄_ell` の構成（スタック・因子・直線束）は本プロジェクトの語彙の外。★無限遠因子は `j` の極因子なので `ht_∞ = h(j)` が標準の読みであり、★★**後続（`Lemma 3.5`・`Lemma 3.7`・`Cor 4.3/4.4`）は 3 本の不等式しか使わない**ので、この読み替えは消費側に影響しない |
| `M_ell(ℚ̄)^{≤d}` の点の集合 | モジュライ点 | ★**`j` 不変量の像** | ℚ̄ 上の楕円曲線は `j` で同型類が決まる。★「点の集合が有限」＝「`j` の集合が有限」 |
| `≲` の向き | 印字は `≲` | ★**`BDge`** | ★`Check/GenEll/Section3NotProvable.lean` の `bdle_chain_forces_bounded` が示すとおり、印字どおり（`BDle`）に読むと `ht_∞` が上に有界になってしまう。`Skeleton` の `prop_3_4` と**同じ読み** |
| semi-stable | 原文の仮定 | ★`SemistableAt`（**受ける**） | ★原文自身の仮定である。★★`primeSubring p` の DVR 付値と `p.valuation L` を繋ぐ橋が無いので、mathlib の `HasMultiplicativeReduction` ではなく本プロジェクトの `valAdd` の言葉で書いている |
| 族の添字 | モジュライ点 | ★添字型 `P` | `Found/GenEll/Prop17.lean` と同じ流儀（点にわたる 1 本の主張） |

★★★★★★**空虚ではない**——3 本の不等式も有限性も、
`§9-1000`〜`§9-1008` で**実際に証明された**ものを使っている。
受けているのは**原文自身の仮定**（半安定・次数の上界）だけである。
-/

namespace ABC3.Found.GaloisRep

open NumberField WeierstrassCurve ABC3.Found.GenEll

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★`Proposition 3.4` -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★**[GenEll] Proposition 3.4**
(Faltings Heights and the Divisor at Infinity)。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

    `deg_∞ ≳ ht_∞ ≳ 12(1+ϵ)·ht^Falt ≳ (1+ϵ)·ht_∞`
    `{ht^Falt ≤ C}` の `j` の像は有限

★仮定は**原文自身のもの**だけである:

| 仮定 | 原文 |
|---|---|
| `hss` | 「with **semi-stable reduction** at all the finite primes」 |
| `hdeg` | `M_ell(ℚ̄)^{≤d}`（定義体の次数が `d` 以下） |

★★中身の出どころ:

| 主張 | 出どころ |
|---|---|
| 第 1 の `≲` | `§9-1008` `degInfOf_le_htFinJ` ＋ `htArchJ_nonneg` |
| 第 2 の `≲` | ★**無条件** `§9-1000`〜`§9-1003` `exists_htJ_le_htFalt'` |
| 第 3 の `≲` | `§9-1006`〜`§9-1007` `exists_htFalt_le_htJ` |
| 有限性 | ★**無条件** `§9-1005` `finite_j_of_htFalt_le` |

★★★逸脱（`ht_∞ ≔ h(j)`、点の集合 ≔ `j` の像、`≲` ≔ `BDge`）は
ファイル冒頭の表に記録した。 -/
theorem prop_3_4 (eps : ℝ) (heps : 0 < eps) (d : ℕ) {P : Type}
    (fld : P → IntermediateField ℚ ℂ) (hnf : ∀ p, NumberField (fld p))
    (hdeg : ∀ p, Module.finrank ℚ (fld p) ≤ d)
    (E : ∀ p, WeierstrassCurve (fld p)) (hE : ∀ p, (E p).IsElliptic)
    (hss : ∀ p, haveI := hnf p; haveI := hE p; ∀ q, SemistableAt q (E p)) :
    BDge (fun p => haveI := hnf p; haveI := hE p; degInfOf (fld p) (E p))
         (fun p => haveI := hnf p; haveI := hE p; htJ (fld p) (E p))
  ∧ BDge (fun p => haveI := hnf p; haveI := hE p; htJ (fld p) (E p))
         (fun p => haveI := hnf p; haveI := hE p; 12 * (1 + eps) * htFaltOf (fld p) (E p))
  ∧ BDge (fun p => haveI := hnf p; haveI := hE p; 12 * (1 + eps) * htFaltOf (fld p) (E p))
         (fun p => haveI := hnf p; haveI := hE p; (1 + eps) * htJ (fld p) (E p))
  ∧ (∀ C : ℝ, ((fun p : P => haveI := hE p; (((E p).j : fld p) : ℂ))
      '' {p | haveI := hnf p; haveI := hE p; htFaltOf (fld p) (E p) ≤ C}).Finite) := by
  obtain ⟨C, hC⟩ := prop_3_4_chain_semistable eps heps
  refine ⟨⟨0, fun p => ?_⟩, ⟨C, fun p => ?_⟩, ⟨C, fun p => ?_⟩, fun C₀ => ?_⟩
  · haveI := hnf p; haveI := hE p
    have h := (hC (fld p) (E p) (hss p)).1
    simp only
    linarith
  · haveI := hnf p; haveI := hE p
    have h := (hC (fld p) (E p) (hss p)).2.1
    simp only
    linarith
  · haveI := hnf p; haveI := hE p
    have h := (hC (fld p) (E p) (hss p)).2.2
    simp only
    linarith
  · exact finite_j_of_htFalt_le d fld hnf hdeg E hE
      (fun p => haveI := hnf p; haveI := hE p; htFaltOf (fld p) (E p)) (fun p => rfl) C₀

/-! ## ★出典の紐付け(`.src`)——★★**項目まるごと** -/

def prop_3_4.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17, item := "Proposition 3.4",
    sectionId := "genell-prop-3-4" }

def prop_3_4.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_htJ_le_htFalt'(第 2 の ≲——無条件、§9-1003)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_htJ_le_htFalt'") 4,
    .citation "[ABC3]" "exists_htFalt_le_htJ(第 3 の ≲、§9-1007)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_htFalt_le_htJ") 4,
    .citation "[ABC3]" "degInfOf_le_htFinJ(第 1 の ≲、§9-1008)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.degInfOf_le_htFinJ") 3,
    .citation "[ABC3]" "finite_j_of_htFalt_le(有限性——無条件、§9-1005)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.finite_j_of_htFalt_le") 4,
    .citation "[mathlib]" "NumberField.logHeight₁_eq(j の Weil 高さの素点分解)"
      (.inMathlib "NumberField.logHeight₁_eq") 2,
    .implicitStep
      ("★★★★★★★★逸脱(2026-08-29): 原文の ht_∞ は M̄_ell の**無限遠因子に付随する高さ**である。" ++
       "本ファイルはそれを **j の素朴 Weil 高さ h(j)** と読み替えた。" ++
       "★M̄_ell の構成(スタック・因子・直線束)は本プロジェクトの語彙の外であり、" ++
       "無限遠因子は j の極因子なのでこれが標準の読みである。" ++
       "★★★**後続(Lemma 3.5・Lemma 3.7・Cor 4.3/4.4)は 3 本の不等式しか使わない**ので、" ++
       "この読み替えは消費側に影響しない") 8,
    .implicitStep
      ("★★★逸脱(2026-08-29): 原文の「M_ell(ℚ̄)^{≤d} の点の集合が有限」を" ++
       "**j 不変量の像が有限**と読んだ。ℚ̄ 上の楕円曲線は j で同型類が決まるので同じ内容である") 6,
    .implicitStep
      ("★★★★逸脱(2026-08-26 に登記済み): 3 本の ≲ を**印字どおり BDle で読むと ht_∞ が" ++
       "上に有界になってしまう**(Check/GenEll/Section3NotProvable.lean の " ++
       "bdle_chain_forces_bounded)。Skeleton/GenEll/Section3.lean と**同じく BDge** で書いた") 7,
    .implicitStep
      ("☆半安定性は原文自身の仮定として受けている(SemistableAt)。" ++
       "★mathlib の HasMultiplicativeReduction ではなく本プロジェクトの valAdd の言葉で" ++
       "書いているのは、primeSubring p の DVR 付値と p.valuation L を繋ぐ橋がまだ無いため" ++
       "(Found/GaloisRep/SemistableFin.lean に記録)") 6,
    .implicitStep
      ("★★★★★★Skeleton からの前進(2026-08-29): Skeleton の prop_3_4 は " ++
       "EllModuliData の 4 つの公理(degInf_le_htInf・htInf_bdeq_faltings・" ++
       "faltingsHeight_bddBelow・northcott)から**組み立てるだけ**だった。" ++
       "★本ファイルはその中身を具体的な量で作る——§9-1000〜§9-1008 の 9 ブロックがそれである。" ++
       "★★定数はすべて ϵ にしか依らない(族にも体にも曲線にも依らない)") 9 ]

end ABC3.Found.GaloisRep
