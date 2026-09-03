/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.VeluSetCard
import ABC3.Found.GenEll.EllModuliObjects
import ABC3.Found.GaloisRep.Lemma35Concrete
import ABC3.Found.GaloisRep.UniformFamily
import ABC3.Found.GaloisRep.BadPrimeData
import ABC3.Meta.Claim

/-!
# QuotClassExists —— `[GenEll] Lemma 3.5` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep IsDedekindDomain NumberField ABC3.Meta
open scoped Classical

/-- ★★★★★★★★★★★★
**商の類は存在する**——★（第 1241）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆位数 `l` の有理点 `Q` と、その Vélu の商が楕円・半安定・乗法還元であることから
`IsQuotClassJ` の存在が出る。★個数の条件は第 1240。

★★逸脱: 商が乗法還元を持つことは仮説として受ける
（原文は「同種なので自動」と括弧で述べる）。 -/
theorem exists_isQuotClassJ (x : RealizedClass) {l : ℕ} (E : SSCurve)
    (hEj : E.j = x.cls) (hEpr : E.PrimeToLocalHeights l)
    (Q : E.W.toAffine.Point) (hQ : addOrderOf Q = l)
    (hell : (veluQuotientFull E.W
      (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))).IsElliptic)
    (hss : ∀ p : HeightOneSpectrum (𝓞 E.fld),
      SemistableAt p (veluQuotientFull E.W
        (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))))
    (hmult : (quotSSCurve E
      (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))
      hell hss).HasMultRed) :
    ∃ y : RealizedClass, IsQuotClassJ x l y.1 := by
  refine ⟨⟨(quotSSCurve E
      (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))) hell hss).j,
    ⟨⟨quotSSCurve E
      (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))) hell hss,
      hmult⟩, rfl⟩⟩, ?_⟩
  exact ⟨E, hEj, hEpr, Q, hQ, hell, hss, rfl⟩

/-! ## ★★★★★★★★★★★★商の類の `deg∞` -/

/-- ★★★★★★★★★★★★
**商の類の `deg∞` は `l` 倍**——★**無条件**（第 1242）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`degInfJ_eq`（類の値は曲線の `deg∞`）と
`degInfOf_eq_of_local`（局所の関係を足し上げる、在庫）を繋いだだけ。

★★★これが `EllModuliWitness` の `degInfJ_quotLCyclicJ` が要る形である
——局所の関係 `v_p(Δ_min(E′)) = l·v_p(Δ_min(E))` は `Lemma 3.2, (ii)`。 -/
theorem degInfJ_quot_eq (E : SSCurve) (S : Finset (E.fld × E.fld))
    (hell : (veluQuotientFull E.W S).IsElliptic)
    (hss : ∀ p : HeightOneSpectrum (𝓞 E.fld), SemistableAt p (veluQuotientFull E.W S))
    (l : ℕ)
    (hloc : ∀ p : HeightOneSpectrum (𝓞 E.fld),
      minDeltaExp p (veluQuotientFull E.W S) = l * minDeltaExp p E.W) :
    degInfJ (quotSSCurve E S hell hss).j = (l : ℝ) * degInfJ E.j := by
  rw [degInfJ_eq, degInfJ_eq]
  exact ABC3.Found.GaloisRep.degInfOf_eq_of_local E.W (veluQuotientFull E.W S) l hloc

/-! ## ★★★★★★★★★★★★商の類の `ht^Falt` -/

/-- ★★★★★★★★★★★★
**商の類の `ht^Falt` の評価は曲線の評価そのもの**——★**無条件**（第 1243）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`faltingsHeightJ_eq`（類の値は曲線の `ht^Falt`、在庫）で置き換えるだけ。

★★★これが `EllModuliWitness` の `faltingsHeightJ_quotLCyclicJ` が要る形である
——同種写像の高さ評価 `hfalt` をそのまま類の言葉に移す。 -/
theorem faltingsHeightJ_quot_le (E : SSCurve) (S : Finset (E.fld × E.fld))
    (hell : (veluQuotientFull E.W S).IsElliptic)
    (hss : ∀ p : HeightOneSpectrum (𝓞 E.fld), SemistableAt p (veluQuotientFull E.W S))
    (l : ℕ) (C₀ : ℝ)
    (hfalt : htFaltOf E.fld (veluQuotientFull E.W S)
      ≤ htFaltOf E.fld E.W + 2 * Real.log l + C₀) :
    faltingsHeightJ (quotSSCurve E S hell hss).j
      ≤ faltingsHeightJ E.j + 2 * Real.log l + C₀ := by
  rw [faltingsHeightJ_eq, faltingsHeightJ_eq]
  exact hfalt

/-! ## ★★★★★★★★★★★★★★★★★★★★第 1339 —— `ht^Falt` は無条件になった -/

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`quotLCyclicJ` の `ht^Falt` の評価**——★**無条件**（第 1339）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆商の類が存在すれば第 1338（無条件の同種写像の高さ評価）がそのまま効き、
存在しなければ `quotLCyclicJ x l = x` なので `2·log l ≥ 0` で済む。

★★★これが `EllModuliWitness` の `faltingsHeight_quotLCyclic` **そのもの**である
——仮説は `l` が素数だけである。 -/
theorem faltingsHeightJ_quotLCyclicJ_uncond (x : RealizedClass) {l : ℕ} (hl : l.Prime) :
    faltingsHeightJ (quotLCyclicJ x l).cls
      ≤ faltingsHeightJ x.cls + 2 * Real.log l := by
  have hlog : (0 : ℝ) ≤ Real.log l :=
    Real.log_nonneg (by exact_mod_cast hl.one_lt.le)
  by_cases hex : ∃ y : RealizedClass, IsQuotClassJ x l y.1
  · obtain ⟨E, hEj, hEpr, Q, hQ, hell, hss, hj⟩ := quotLCyclicJ_spec x l hex
    haveI := hell
    have hfalt := ABC3.Found.GaloisRep.htFalt_veluQuotientFull_le_uncond
      E.W
      (veluQuotientFull E.W
        (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))))
      hl Q hQ rfl
    have hq := faltingsHeightJ_quot_le E _ hell hss l 0 (by linarith [hfalt])
    rw [hj, hEj] at hq
    show faltingsHeightJ (quotLCyclicJ x l).1 ≤ _
    linarith [hq]
  · have hself : quotLCyclicJ x l = x := by
      rw [quotLCyclicJ, dif_neg hex]
    rw [hself]
    linarith

/-! ## ★出典の紐付け(`.src`)——★**条つきである。指標には数えない** -/

def faltingsHeightJ_quotLCyclicJ_uncond.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(quotLCyclicJ の ht^Falt の評価。★無条件)",
    sectionId := "genell-lemma-3-5" }

def faltingsHeightJ_quotLCyclicJ_uncond.needs : List ProofObligation :=
  [ .citation "[ABC3]" "htFalt_veluQuotientFull_le_uncond(第 1338、無条件)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.htFalt_veluQuotientFull_le_uncond") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1339）**——`EllModuliWitness` の " ++
       "`faltingsHeight_quotLCyclic` **そのもの**である。" ++
       "☆`IsQuotClassJ` を生成元 `Q` を持ち歩く形に締め直したのが効いた。") 3 ]

def faltingsHeightJ_quotLCyclicJ_of_isog.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(quotLCyclicJ の ht^Falt の評価——同種写像の高さ評価を仮説で受ける)",
    sectionId := "genell-lemma-3-5" }

def degInfJ_quotLCyclicJ_of_local.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(quotLCyclicJ の deg∞ は l 倍——局所の関係は Lemma 3.2, (ii) を仮説で受ける)",
    sectionId := "genell-lemma-3-5" }

def faltingsHeightJ_quot_le.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(商の類の ht^Falt の評価は曲線の評価そのもの。★無条件)",
    sectionId := "genell-lemma-3-5" }

def degInfJ_quot_eq.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(商の類の deg∞ は l 倍。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_isQuotClassJ.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(商の類の存在——商の乗法還元は仮説として受ける)",
    sectionId := "genell-lemma-3-5" }

def exists_isQuotClassJ.needs : List ProofObligation :=
  [ .citation "[ABC3]" "card_image_pointCoords_nsmul(第 1240、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.card_image_pointCoords_nsmul") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1241）**——`quotLCyclicJ_spec`（`EllModuliObjects`）を" ++
       "使えるようにする段である。☆★★逸脱: 商が**乗法還元を持つ**ことは" ++
       "仮説として受ける（原文は「同種なので自動」と括弧で述べる。" ++
       "`VeluQuotOK` と同じ扱い）。") 3 ]

end ABC3.Found.GenEll
