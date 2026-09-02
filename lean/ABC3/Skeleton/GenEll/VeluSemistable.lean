/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.VeluQuotElliptic
import ABC3.Meta.Claim

/-!
# 第 1342 ブロック —— **Vélu の商の半安定性**（`Skeleton`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★★★★★これは何か——**残った 1 本を節点にする**

原文が「同種なので自動」と括弧で述べる段は 2 枚に分かれていた:

| 枚 | 内容 | 状態 |
|---|---|---|
| 楕円性 | `veluQuotientFull E ⟨Q⟩` が楕円 | ★**閉じた**（第 1330-1336） |
| 半安定性 | 各素点で半安定 | ☆悪い素点は第 1327、**良い素点が残る** |

★★★本ファイルは**残った半安定性を 1 つの節点として立てる**。
これで `VeluQuotOK`（したがって `Lemma 3.7` の安定直線の側と
`EllModuliWitness` の `lcyclicExc` 欄）が**この 1 本だけ**に依る形になる。

## ★★★★★★★★何が要るのか——測定（2026-09-02、第 1342）

☆**悪い素点**（`v_p(j) < 0`）——`semistableAt_velu_of_veluCurve_eq`（第 1327、証明済み）。
残るのは Tate モデルへの局所の配管である。

★**良い素点**（`v_p(j) ≥ 0`）——`minDeltaExp p E′ = 0`、すなわち
**同種で良還元が保たれること**が要る。道は

1. `E` の極小モデルは `v_p(Δ) = 0`
2. `l`-捩れ点の座標は整（`mem_primeSubring_x_of_addOrderOf_prime`、第 1073、証明済み）
3. Vélu の `v`・`w` は整なので商も整（`isIntegral_veluQuotientFull_of_addOrderOf_prime`、第 1074）
4. **剰余体の上で `Δ(Ẽ/H̃) ≠ 0`**

であり、根は **(4) の「剰余体の上の Vélu の定理」**である。

★★☆**なぜ `ℂ` の道が使えないか**——第 1330-1335 は一意化（`℘` 函数）で
Vélu の定理を取った。☆これは標数 `0` の解析的な道であり、標数 `p` の剰余体では使えない。
★別の道（Néron–Ogg–Shafarevich、あるいは剰余体上の代数的な Vélu の定理）が要る。
-/

namespace ABC3.Skeleton.GenEll

open WeierstrassCurve IsDedekindDomain NumberField Finset
open ABC3.Found.GenEll ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] Vélu の商は半安定**（第 1342、第 1345 で一般形に直した）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆原文が「同種なので自動」と括弧で述べる段の**半安定性の側**である。

★★★これが `VeluQuotOK`・`IsQuotClassJ` に残る**ただ 1 つの節点**である。 -/
theorem semistableAt_veluQuotientFull {L : Type} [Field L] [NumberField L] [DecidableEq L]
    (E : WeierstrassCurve L) [E.IsElliptic]
    (hss : ∀ p : HeightOneSpectrum (𝓞 L), SemistableAt p E)
    {l : ℕ} (Q : E.toAffine.Point) (hQ : addOrderOf Q = l)
    (p : HeightOneSpectrum (𝓞 L)) :
    SemistableAt p (veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) := by
  sorry

/-- ★★★★★★★★★★★★**`SSCurve` の語彙で**（第 1345）。 -/
theorem semistableAt_veluQuot_ss (E : SSCurve) {l : ℕ} (Q : E.W.toAffine.Point)
    (hQ : addOrderOf Q = l) (p : HeightOneSpectrum (𝓞 E.fld)) :
    SemistableAt p (veluQuotientFull E.W
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) :=
  semistableAt_veluQuotientFull E.W E.ss Q hQ p

/-- ★★★★★★★★★★★★★★★★★★★★
**`VeluQuotOK` はこの 1 本から出る**——★（第 1342）。

☆楕円性は第 1336（無条件）、半安定性は上の節点。 -/
theorem veluQuotOK_all (E : SSCurve) (l : ℕ) : VeluQuotOK E l := by
  refine veluQuotOK_of_semistable E l (fun M _ Q' hQ' P => ?_)
  letI : DecidableEq (M : Type) := fun a b => Classical.propDecidable (a = b)
  letI : NumberField (M : Type) := NumberField.of_module_finite E.fld M
  haveI hEell : (E.W.baseChange (M : Type)).IsElliptic := by
    show (E.W.map (algebraMap E.fld M)).IsElliptic
    infer_instance
  exact semistableAt_veluQuotientFull (E.W.baseChange (M : Type))
    (ABC3.Found.GaloisRep.semistableAt_baseChange_all E.fld (M : Type) E.W E.ss) Q' hQ' P

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def semistableAt_veluQuotientFull.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の商は半安定——原文が「同種なので自動」と括弧で述べる段)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_veluQuotientFull.needs : List ProofObligation :=
  [ .citation "[ABC3]" "semistableAt_velu_of_veluCurve_eq(悪い素点側、第 1327、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.semistableAt_velu_of_veluCurve_eq") 1,
    .citation "[ABC3]" "isIntegral_veluQuotientFull_of_addOrderOf_prime(良い素点で商は整、第 1074)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.isIntegral_veluQuotientFull_of_addOrderOf_prime") 1,
    .citation "[FC]" "Degenerations of Abelian Varieties(同種で還元型が保たれること)"
      (.absent
        ("mathlib に Néron モデル・導手・Néron–Ogg–Shafarevich は無い(2026-09-02 確認)。" ++
         "★本プロジェクトが `ℂ` で取った Vélu の定理(第 1330-1335)は一意化を使うので" ++
         "標数 p の剰余体では使えない——剰余体上の代数的な Vélu の定理が要る")) 17,
    .implicitStep
      ("★★★★**2026-09-02（第 1342）**——原文は括弧 1 つで済ませている段である。" ++
       "☆楕円性の側は第 1330-1336 で閉じた。悪い素点の半安定性は第 1327 で閉じた。" ++
       "★残るのは**良い素点で `minDeltaExp p E′ = 0`**（同種で良還元が保たれること）" ++
       "ただ 1 つであり、その根は剰余体の上の Vélu の定理 `Δ(Ẽ/H̃) ≠ 0` である。") 17 ]

def veluQuotOK_all.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(VeluQuotOK はこの 1 本から出る——楕円性は第 1336 で無条件)",
    sectionId := "genell-lemma-3-5" }

def veluQuotOK_all.needs : List ProofObligation :=
  [ .citation "[ABC3]" "veluQuotOK_of_semistable(第 1336、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.veluQuotOK_of_semistable") 1,
    .citation "[ABC3]" "semistableAt_veluQuotientFull(本ファイルの節点)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.semistableAt_veluQuotientFull") 1 ]

def semistableAt_veluQuot_ss.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の商は半安定——SSCurve の語彙で)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_veluQuot_ss.needs : List ProofObligation :=
  [ .citation "[ABC3]" "semistableAt_veluQuotientFull(本ファイルの節点)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.semistableAt_veluQuotientFull") 1 ]

end ABC3.Skeleton.GenEll
