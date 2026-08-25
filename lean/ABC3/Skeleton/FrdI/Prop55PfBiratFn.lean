import ABC3.Meta.Claim
import ABC3.Found.FrdI.Prop55ScaleRootCoa
import ABC3.Found.FrdI.Prop55Std

/-!
# [FrdI] Proposition 5.5, (iii) —— `𝒞^pf` は birationally Frobenius-normalized 型(`Skeleton`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.105。

原文 (FrdI p.105):
> if C is of standard (respectively, rationally standard) type, then so are Cun-tr, Crlf.

## ★★このファイルの位置づけ —— `Theorem 6.4, (i)` の律速 1 点

`Theorem 6.4, (i)` は `C, C^pf, C^rlf, C^un-tr, (C^pf)^un-tr` の 5 圏が
rationally standard だと言う。★2026-08-25 現在:

| 圏 | 状態 |
|---|---|
| `𝒞` | 済(`ex61_ratStd` / `ex63_ratStd`) |
| `𝒞^un-tr` | 済(`unTr_isOfRationallyStandardType`) |
| `𝒞^pf` | ★**本ファイルの 1 条だけが残っている** |
| `𝒞^rlf` | 未 |
| `(𝒞^pf)^un-tr` | 未(`(𝒞^un-tr)^pf` 経由の道がある) |

★`𝒞^pf` の `Definition 4.5, (iii)` の 4 条のうち
**standard**(`pfRoot_standardType`)・
**`((𝒞^pf)^un-tr)^birat` の Frobenius-compact 対象**(`unTr_pf_biratCompact_of_baseId`)・
**rational の台の対応**(`Def24PfTransport.lean`)は済んでおり、
残るのが**本ファイルの birationally Frobenius-normalized** である。

## ★★★なぜ在庫の `pfRoot_isOfModelType` では出ないか

在庫の `pfRoot_isOfModelType_of_arb` は 4 条を要求し、そのうち
**unit-trivial 型**だけが `Example 6.1` / `Example 6.3` の `𝒞` で成り立たない
(単元が `k_L^×` / `μ(L)`)。★これが `hut` の壁である。

## ★★★★手順書(2026-08-25 に測った。次に着手する人へ)

★**`isFrobeniusNormalized_transport`(在庫)は関手を要求しない** ——
1 つの対象での `End` の乗法同型 ＋ 3 つの対応だけである。

1. **`End` の全単射は在庫**。`biratPfHomEquivRoot`(`Prop55ScaleRootCoa.lean`)を
   `A = B`、`n = m`、`k = n·n` で読むと

       End_{(𝒞^birat)^pf}(biratUp (rtObj A n))  ≃  End_{(𝒞^pf)^birat} ⟨A, n⟩

   がそのまま出る(★2026-08-25 に型検査で確認済み)。
2. **乗法性 —— ★2026-08-25 に閉じた**(`Found/FrdI/Prop55ScaleRootBirat.lean`):

       biratPfEndMulEquiv        (1 段目、`biratPfHom_comp` が乗法性)
       scaleRootBiratEndMulEquiv (2 段目、`Functor.mapEnd` は `MonoidHom`)
       pfBiratEndMulEquiv        (3 段を束ねた `End` の乗法同型)

   ★3 段目の共役は `Iso.conj`(mathlib、`End X ≃* End Y`)を使ったので乗法性は無料。
3. **3 つの対応**(`IsBaseIdentity` / `degFr` / `OTri`):
   * 2 段目(`Σ_k` の birat 版)—— ★**2026-08-25 に閉じた**
     (`biratBase_scaleRootBirat` / `biratDeg_scaleRootBirat`)。
   * 1 段目(`biratPfHom`)—— 在庫の
     `biratBase_biratPfHom` / `biratDeg_biratPfHom` / `otimes_biratPfHom` は
     **引数が `toHomPf ε` の形のものだけ**である。
     ★`End` の一般の元は `HomPf.mk W z`(`homPf_birat_exists_rep`)なので、
     **`pfBase_mk` / `pfDeg_mk`(`repBase` / `repDeg`)経由の一般版が要る** ——
     ★★**ここが本ファイルに残る唯一の実作業**である。
   * 3 段目(同型による共役)—— 一般論で済む。
4. 仕上げは `pf_frobNormalizedType`(`𝒞^birat` 側、在庫)を流すだけ。
   ★`IsOfFrobeniusNormalizedType (biratPre P G)` は
   `IsOfBirationallyFrobeniusNormalizedType C P G` と**同じもの**である
   (`BiratCat P G = C` かつ `(toBiratCat).obj A = A` なので)。

★★見積りの推移(2026-08-25 の 1 日で 3 回下がった):
**400 行 → 150 行 → 上の (2)(3) だけ**。
教訓は毎回同じ: **在庫を検索してから見積もる。**
-/

namespace ABC3.Skeleton.Prop55

open CategoryTheory ABC3.Found.FrdI ABC3.Meta

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] [IsConnected D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}
  {G : Frobenioid P}

/-- ★★★★★★**[FrdI] Proposition 5.5, (iii)** ——
`𝒞` が birationally Frobenius-normalized 型なら **`𝒞^pf` もそう**。

原文 (FrdI p.105):
> if C is of standard (respectively, rationally standard) type, then so are Cun-tr, Crlf.

★★これが `Theorem 6.4, (i)` の 5 圏のうち `𝒞^pf` を止めている**唯一の条**である。
手順書はファイル冒頭の「★★★★手順書」を見よ。

★仮定は原文の常備仮定(`hfi` / `hiso`)と `𝒞` の birat-Frobenius-normalized 性だけで、
**`hut`(unit-trivial)は要らない**——そこが在庫の `pfRoot_isOfModelType_of_arb` との違いである。 -/
theorem pfRoot_biratFrobNormalizedType
    (hfi : IsOfFrobeniusIsotropicType P) (hiso : ∀ X : C, IsIsotropic P X)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G))
    (hbfn : IsOfBirationallyFrobeniusNormalizedType C P G) :
    IsOfBirationallyFrobeniusNormalizedType (PfRootObj P F) (pfRootPre P F) Gpf := by
  sorry

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def pfRoot_biratFrobNormalizedType.src : Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — 𝒞^pf は birationally Frobenius-normalized 型",
    sectionId := "frdi-prop-5-5" }

def pfRoot_biratFrobNormalizedType.needs : List ProofObligation :=
  [ .citation "[ABC3]" "biratPfHomEquivRoot(一般の根での Hom の全単射。A=B・n=m・k=n·n で End になる)"
      (.inProject "ABC3" "ABC3.Found.FrdI.biratPfHomEquivRoot") 105,
    .citation "[ABC3]" "isFrobeniusNormalized_transport(End の乗法同型 ＋ 3 つの対応だけで移る)"
      (.inProject "ABC3" "ABC3.Found.FrdI.isFrobeniusNormalized_transport") 105,
    .citation "[ABC3]" "pf_frobNormalizedType(𝒞^birat の側。これが移送の出発点)"
      (.inProject "ABC3" "ABC3.Found.FrdI.pf_frobNormalizedType") 105,
    .citation "[ABC3]" "biratPfHom_id / biratPfHom_comp(1 段目の乗法性)"
      (.inProject "ABC3" "ABC3.Found.FrdI.biratPfHom_comp") 105,
    .citation "[ABC3]" "biratBase_biratPfHom / biratDeg_biratPfHom / otimes_biratPfHom(1 段目の 3 対応)"
      (.inProject "ABC3" "ABC3.Found.FrdI.biratBase_biratPfHom") 105,
    .citation "[ABC3]" "scaleRootBirat(2 段目。根の一斉倍化の birat 版、関手なので乗法的)"
      (.inProject "ABC3" "ABC3.Found.FrdI.scaleRootBirat") 105,
    .citation "[ABC3]" "rootBase_scaleRootHom / rootDeg_scaleRootHom(2 段目の 3 対応の材料)"
      (.inProject "ABC3" "ABC3.Found.FrdI.rootBase_scaleRootHom") 105,
    .citation "[ABC3]" "pfRoot_exists_iso_root(3 段目の共役に使う同型 ⟨A,n⟩ ≅ ⟨A^{(m)}, n·m⟩)"
      (.inProject "ABC3" "ABC3.Found.FrdI.pfRoot_exists_iso_root") 105,
    .implicitStep
      ("★実作業は 2 つだけ: (2) 3 段の合成が乗法的であること、" ++
       "(3) 2 段目(scaleRootBirat)の Base・degFr・OTri の保存を " ++
       "biratPsiMap_mk で代表元に降ろして示すこと。" ++
       "★1 段目の 3 対応と 3 段目の共役は在庫と一般論で済む") 105 ]

end ABC3.Skeleton.Prop55
