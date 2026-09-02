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

## ★★★★★★★★★★★★★★★★道の候補を測った（2026-09-02、第 1384）

★★★**これが `EllModuliWitness` に残るただ 1 つの葉になった**
——α 側（`exists_h2_h1_unipotent`）は第 1384 で閉じた。

☆良い素点で既にあるもの:

* `E` が `p` 上整で `p ∤ l` なら **Vélu の商も `p` 上整**（第 1074、証明済み）
* `Δ(E′) ≠ 0`（大域の楼円性、第 1336、無条件）

★残るのは **`v_p(Δ(E′)) = 0`** だけである。

### ☆道の候補 A——**普遍な整除性 `Δ(E′) ∣ Δ(E)^l`**

★`l = 2` で手計算して確かめた（2026-09-02）——
`E : y² = x³ + ax + b`、核 `{O, (x₀,0)}`、`t = 3x₀² + a` のとき

    Δ(E′) = Δ(E)² / t³        すなわち  Δ(E)² = Δ(E′)·t³

（検算: `a = −1, b = 0, x₀ = 0` なら `t = −1`、`E′ : y² = x³+4x`、
`Δ = 64`、`Δ′ = −4096 = 64²/(−1)³`。）

★★**この形の普遍恒等式 `Δ(E)^l = Δ(E′)·G`（`G` は核の座標の多項式）
があれば、良い素点では `v_p(Δ(E)^l) = 0` から
直ちに `v_p(Δ(E′)) = 0` が出る**（両因子が整なので）。
☆`G` の**明示の形は要らない**のがこの道の利点である。

☆さらに、恒等式は**係数の多項式恒等式**なので、
本プロジェクトが既に使っている「一般の媒介変数で書いて `ℂ` で確かめる」
（`genericParam`・`c4DenomFree` の手）が使える可能性がある。

### ☆道の候補 B——剰余体の上の代数的な Vélu

★「核が位数 `l` の部分群なら `Δ(E/C) ≠ 0`」を任意の体で直接取る道。
☆商の滑らかさ（有限部分群スキームによる商）が要り、mathlib には無い。

### ☆道の候補 C——Néron–Ogg–Shafarevich

★mathlib に Néron モデル・導体・ NOS は無い（2026-09-02 確認）。

★★★**候補 A が最も短い**——他の 2 つは新しい理論を丸ごと建てることになる。

### ★★★★★★★★候補 A の**恒等式を見つけた**（2026-09-02、第 1384）

    Δ(E)^l = Δ(E/C) · ( ∏_{P ∈ C∖{O}} (2 y_P + a₁ x_P + a₃) )^4

☆`l = 3, 5, 7` について複数の曲線（Tate 標準形の族）で**厳密に成立**することを
確かめた——`tools/velu-disc-check.py`。

★重みの検算: `Δ` は重み 12、`2y_P + a₁x_P + a₃` は重み 3 なので
右辺の第 2 因子は `3(l−1)·4 = 12(l−1)`、合計 `12 + 12(l−1) = 12l` ✓。

★★★**この恒等式があれば良い素点の葉は即座に閉じる**——
`v_p(Δ(E)) = 0` で、`Δ(E′)` も `∏(2y_P+a₁x_P+a₃)` も整（第 1073-1074）なので、
積が単元 ⇒ **両方とも単元** ⇒ `v_p(Δ(E′)) = 0` ✓。

☆**証明の道**——両辺は体の元なので、
`ℂ` で証明して埋め込みで降ろす手（第 1334 の `isElliptic_velu_congr_curve` と同じ型）が使える。
★`ℂ` では一意化（第 1330-1335）があるので、
`Δ(Λ)` の積公式と `℘′(u) = −σ(2u)/σ(u)⁴` から出るはずである。
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
